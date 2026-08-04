#!/usr/bin/env python3
"""Chromosome-parallel PCR-duplicate removal with pluggable dedup backends.

Rationale (see: markdown docs in repo):
  * UMICollapse deduplicates independently per alignment position.
  * Read pairs always map to the same chromosome.
  * => each chromosome can be deduplicated fully independently, then the
    per-chromosome results merged with `samtools merge`.

This gives real multi-process / wall-clock parallelism while keeping each
per-chromosome dedup job's memory small and bounded.

Flow:
  1. Read @SQ contig list from the input BAM header.
  2. Slice the BAM into one small BAM per contig (streamed via `samtools view`).
  3. Run the selected dedup backend on each slice in parallel.
  4. Merge per-chromosome BAMs -> coordinate-sorted output (+ .bai).
  5. Optional JSON report: per-contig timing / return code + totals.

Backends (extensible via the method table):
  umicollapse-bktree    UMICollapse + BK-tree + two-pass      (recommended)
  umicollapse-fenwick   UMICollapse + Fenwick BK-tree + twopass
  umicollapse-naive     UMICollapse + naive (original behaviour, needs huge heap)
  picard                Picard MarkDuplicates (RNA-safe: no optical detection)

Example:
  dedup_split.py -i combined.bam -o deduped.bam -m umicollapse-bktree \\
                 --threads 32 --report dedup_report.json
"""

import argparse
import concurrent.futures
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time


def log(msg):
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()


# --------------------------------------------------------------------------- #
#  Command construction per backend
# --------------------------------------------------------------------------- #
def _cmd_umicollapse(p, data, paired, threads, java_opts):
    cmd = ["java"] + java_opts.split() + [
        "-jar", p["umicollapse_jar"], "bam",
        "-T", str(threads),
        "--remove-chimeric",
        "--data", data,
        "--merge", "avgqual",
        "--two-pass",
    ]
    if paired:
        cmd += ["--paired", "--remove-unpaired"]
    return cmd


def _cmd_picard(p, threads, java_opts):
    return ["java"] + java_opts.split() + [
        "-jar", p["picard_jar"], "MarkDuplicates",
        "--READ_NAME_REGEX", "null",
        "--DUPLICATE_SCORING_STRATEGY", "SUM_OF_BASE_QUALITIES",
        "--REMOVE_DUPLICATES", "true",
        "--VALIDATION_STRINGENCY", "SILENT",
        "--TMP_DIR", p["tmpdir"],
    ]


def build_command(method, p, in_bam, out_bam, threads, java_opts):
    if method.startswith("umicollapse"):
        data = {
            "umicollapse-bktree": "bktree",
            "umicollapse-fenwick": "fenwickbktree",
            "umicollapse-naive": "naive",
        }[method]
        base = _cmd_umicollapse(p, data, p["paired"], threads, java_opts)
        return base + ["-i", in_bam, "-o", out_bam]
    if method == "picard":
        base = _cmd_picard(p, threads, java_opts)
        return base + ["--INPUT", in_bam, "--OUTPUT", out_bam,
                       "--METRICS_FILE", out_bam + ".picard_metrics.txt"]
    raise ValueError(f"Unknown method: {method}")


# --------------------------------------------------------------------------- #
#  Slicing / merging helpers
# --------------------------------------------------------------------------- #
def get_contigs(samtools, bam):
    header = subprocess.run([samtools, "view", "-H", bam],
                            capture_output=True, text=True, check=True).stdout
    contigs = []
    for ln in header.splitlines():
        if not ln.startswith("@SQ"):
            continue
        for f in ln.split("\t"):
            if f.startswith("SN:"):
                contigs.append(f[3:])
    return contigs


def slice_contig(samtools, in_bam, contig, out_bam):
    subprocess.run([samtools, "view", "-b", "-@", "1", in_bam, contig, "-o", out_bam],
                   check=True)


def _slice_one(args):
    """Worker for parallel slicing: slice one contig + index it."""
    samtools, in_bam, contig, sb = args
    slice_contig(samtools, in_bam, contig, sb)
    subprocess.run([samtools, "index", sb], check=True)
    return sb


# --------------------------------------------------------------------------- #
#  Worker (dedup)
# --------------------------------------------------------------------------- #
def _run_one(task):
    cmd, out_bam, contig = task
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    return {
        "contig": contig,
        "ok": proc.returncode == 0,
        "returncode": proc.returncode,
        "elapsed_s": round(time.time() - t0, 2),
        "stderr_tail": (proc.stderr or "")[-1500:],
    }


# --------------------------------------------------------------------------- #
#  Main orchestration
# --------------------------------------------------------------------------- #
def deduplicate_split(in_bam, out_bam, method, p, threads=16, java_opts="-Xms8g -Xmx64g"):
    samtools = p["samtools"]
    if not os.path.exists(in_bam + ".bai") and not os.path.exists(in_bam.replace(".bam", ".bai")):
        log(f"[dedup_split] building index for input: {in_bam}")
        subprocess.run([samtools, "index", in_bam], check=True)
    with tempfile.TemporaryDirectory(prefix="dedsplit_", dir=p["tmpdir"]) as td:
        contigs = get_contigs(samtools, in_bam)
        log(f"[dedup_split] {len(contigs)} contigs; method={method}; threads={threads}")

        # 1. slice (parallel — each contig is I/O-bound, threads work well)
        slice_tasks = [
            (samtools, in_bam, c, os.path.join(td, f"c{idx:02d}.bam"))
            for idx, c in enumerate(contigs)
        ]
        slice_bams = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=p["parallel"]) as ex:
            futs = {ex.submit(_slice_one, t): t for t in slice_tasks}
            for fut in concurrent.futures.as_completed(futs):
                slice_bams.append(fut.result())
        # Preserve contig order so merge ordering is deterministic
        slice_bams.sort()

        # 2. parallel dedup
        tasks = []
        for sb in slice_bams:
            ob = sb + ".dedup.bam"
            cmd = build_command(method, p, sb, ob, threads, java_opts)
            tasks.append((cmd, ob, os.path.basename(sb)))

        results = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=p["parallel"]) as ex:
            futs = {ex.submit(_run_one, t): t for t in tasks}
            for fut in concurrent.futures.as_completed(futs):
                results.append(fut.result())

        failed = [r for r in results if not r["ok"]]
        for r in failed:
            log(f"  FAILED {r['contig']} rc={r['returncode']}:\n{r['stderr_tail']}")
        if failed:
            raise RuntimeError(f"{len(failed)}/{len(contigs)} contigs failed")

        # 3. merge + index
        deduped = [sb + ".dedup.bam" for sb in slice_bams]
        subprocess.run([samtools, "merge", "-f", out_bam] + deduped, check=True)
        subprocess.run([samtools, "index", out_bam, out_bam + ".bai"], check=True)

        return results, out_bam


def main():
    ap = argparse.ArgumentParser(description="Chromosome-parallel dedup (pluggable backends).")
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("-m", "--method", default="umicollapse-bktree",
                    choices=["umicollapse-bktree", "umicollapse-naive",
                             "umicollapse-fenwick", "picard"])
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--java-opts", default="-Xms8g -Xmx64g")
    ap.add_argument("--pair-count", type=int, default=2, help="1=SE 2=PE")
    ap.add_argument("--report", default=None)
    ap.add_argument("--samtools", default=shutil.which("samtools") or "samtools")
    ap.add_argument("--umicollapse-jar", default=os.path.expanduser("~/tools/UMICollapse/umicollapse.jar"))
    ap.add_argument("--picard-jar", default=os.path.expanduser("~/tools/picard.jar"))
    ap.add_argument("--tmpdir", default=None,
                    help="Scratch dir for per-chromosome slices & JVM temp. "
                         "Defaults to $SLURM_TMPDIR, else $TMPDIR, else the output dir.")
    args = ap.parse_args()

    if not args.java_opts.startswith("-server"):
        args.java_opts = "-server " + args.java_opts

    # Pick a fast local scratch dir for intermediates. SLURM nodes have node-local
    # fast disk at $SLURM_TMPDIR -- prefer it; fall back to $TMPDIR, else a cache
    # subdir next to the output (i.e. on the shared filesystem).
    tmpdir = args.tmpdir or os.environ.get("SLURM_TMPDIR") or os.environ.get("TMPDIR")
    if not tmpdir:
        out_dir = os.path.dirname(os.path.abspath(args.output))
        tmpdir = os.path.join(out_dir, ".dedup_cache")
    os.makedirs(tmpdir, exist_ok=True)

    p = {
        "paired": args.pair_count == 2,
        "parallel": max(1, args.threads),
        "samtools": args.samtools,
        "umicollapse_jar": args.umicollapse_jar,
        "picard_jar": args.picard_jar,
        "tmpdir": tmpdir,
    }
    log(f"[dedup_split] scratch/tmp dir: {tmpdir}")

    t0 = time.time()
    results, out = deduplicate_split(
        args.input, args.output, args.method, p,
        threads=args.threads, java_opts=args.java_opts,
    )
    total = time.time() - t0

    if args.report:
        with open(args.report, "w") as f:
            json.dump({"method": args.method, "output": out,
                       "total_time_s": round(total, 2), "contigs": results}, f, indent=2)
    log(f"[dedup_split] DONE method={args.method} output={out} wall={total:.1f}s")


if __name__ == "__main__":
    main()
