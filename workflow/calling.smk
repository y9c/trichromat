# calling
TEMPDIR = Path(
    os.path.relpath(
        config.get("tempdir", os.path.join(workflow.basedir, ".tmp")), workflow.basedir
    )
)
INTERNALDIR = Path("internal_files")

PATH = config["path"]

LIBTYPE = config["libtype"]
WITH_UMI = config.get(
    "with_umi",
    LIBTYPE in ["INLINE", "ECLIP6", "ECLIP10", "TAKARAV3", "SACSEQ", "SACSEQV3"],
)
MARKDUP = config.get("markdup", True)
# force markdup to be true for UMI
# if WITH_UMI:
#     MARKDUP = True
STRANDNESS = config.get("strandness", True)
GENE_NORC = config.get("gene_norc", True)
BASE_CHANGE = config.get("base_change", "A,G")
SPLICE_GENOME = config.get("splice_genome", True)
SPLICE_CONTAM = config.get("splice_contamination", False)
REF = config.get("_REF", {})
READS = config.get("_READS", {})

MOTIF_FLANKING = config.get("motif_flanking", 2)
DROP_CLUSTERED = config.get("drop_clustered", True)


rule build_internal_fa_index:
    input:
        INTERNALDIR / "reference_file/{reftype}.fa",
    output:
        INTERNALDIR / "reference_file/{reftype}.fa.fai",
    threads: 8
    shell:
        """
        {PATH[samtools]} faidx -@ {threads} {input} --fai-idx {output}
        """


rule build_genome_fa_index:
    input:
        REF["genome"][0],
    output:
        REF["genome"][0] + ".fai",
    threads: 8
    shell:
        """
        {PATH[samtools]} faidx -@ {threads} {input} --fai-idx {output}
        """


# hisat:
# ref     pos     strand  convertedBaseCount      unconvertedBaseCount
# pbr:
# #chrom  pos0    ref_base        depth   a       c       g       t       n
# fwd: $1, $2 + 1, +, $7, $5
# rev: $1, $2 + 1, -, $6, $8

# {PATH[samtools]} view -@ {threads} -e "rlen < 100000 && [XM] * 20 <= (qlen-sclen) && [Zf] <= 3 && 3 * [Zf] <= [Zf] + [Yf]" -h {input} | \
#     {PATH[hisat3ntable]} -p {threads} --alignments - --ref {params.fa} --output-name /dev/stdout --base-change {params.basechange} | \
#     cut -f 1,2,3,5,7 | \
#     gzip > {output}


rule pileup_by_fwd_strand:
    input:
        bam=INTERNALDIR / "aligned_bam/{sample}.{reftype}.bam",
        bai=INTERNALDIR / "aligned_bam/{sample}.{reftype}.bam.bai",
        fa=lambda wildcards: (
            REF["genome"][0]
            if wildcards.reftype == "genome"
            else INTERNALDIR / f"reference_file/{wildcards.reftype}.fa"
        ),
        fai=lambda wildcards: (
            REF["genome"][0] + ".fai"
            if wildcards.reftype == "genome"
            else INTERNALDIR / f"reference_file/{wildcards.reftype}.fa.fai"
        ),
    output:
        temp(TEMPDIR / "pileup_by_strand/{sample}.{reftype}.fwd.tsv.gz"),
    params:
        ref_base=BASE_CHANGE.split(",")[0],
        ref_col={"A": 5, "C": 6, "G": 7, "T": 8}[BASE_CHANGE.split(",")[0]],
        alt_col={"A": 5, "C": 6, "G": 7, "T": 8}[BASE_CHANGE.split(",")[-1]],
        motif_flanking=MOTIF_FLANKING,
        motif_center=MOTIF_FLANKING + 1,
        drop_clustered=(
            'and read:tag("Zf")<=3 and 2 * read:tag("Zf") <= read:tag("Yf")'
            if DROP_CLUSTERED
            else ""
        ),
    threads: 16
    shell:
        """
        {PATH[pb]} -t {threads} --mate-fix --max-depth 500000 -k {params.motif_flanking} --pile-expression 'return string.upper(string.sub(pile.ref_base,{params.motif_center},{params.motif_center}))=="{params.ref_base}"' --fasta {input.fa} -e 'return read.bq>20 and read.strand==1 and read:tag("XM") * 20 <= read.length and {params.drop_clustered}' {input.bam} |\
            awk -v OFS="\\t" 'NR>1 && ${params.alt_col}+${params.ref_col}>0 {{print $1,$2+1,"+",$3,${params.alt_col},${params.ref_col}}}' | \
            gzip > {output}
        """


rule pileup_by_rev_strand:
    input:
        bam=INTERNALDIR / "aligned_bam/{sample}.{reftype}.bam",
        bai=INTERNALDIR / "aligned_bam/{sample}.{reftype}.bam.bai",
        fa=lambda wildcards: (
            INTERNALDIR / "reference_file/genes.fa"
            if wildcards.reftype == "genes"
            else REF["genome"][0]
        ),
        fai=lambda wildcards: (
            INTERNALDIR / "reference_file/genes.fa.fai"
            if wildcards.reftype == "genes"
            else REF["genome"][0] + ".fai"
        ),
    output:
        temp(TEMPDIR / "pileup_by_strand/{sample}.{reftype}.rev.tsv.gz"),
    params:
        ref_base={"A": "T", "C": "G", "G": "C", "T": "A"}[BASE_CHANGE.split(",")[0]],
        # a c g t
        # 5,6,7,8
        # into
        # 8,7,6,5
        ref_col={"A": 8, "C": 7, "G": 6, "T": 5}[BASE_CHANGE.split(",")[0]],
        alt_col={"A": 8, "C": 7, "G": 6, "T": 5}[BASE_CHANGE.split(",")[-1]],
        motif_flanking=MOTIF_FLANKING,
        motif_center=MOTIF_FLANKING + 1,
        drop_clustered=(
            'and read:tag("Zf")<=3 and 2 * read:tag("Zf") <= read:tag("Yf")'
            if DROP_CLUSTERED
            else ""
        ),
    threads: 16
    shell:
        """
        {PATH[pb]} -t {threads} --mate-fix --max-depth 500000 -k {params.motif_flanking} --pile-expression 'return string.upper(string.sub(pile.ref_base,{params.motif_center},{params.motif_center}))=="{params.ref_base}"' --fasta {input.fa} -e 'return read.bq>20 and read.strand==-1 and read:tag("XM") * 20 <= read.length {params.drop_clustered}' {input.bam} |\
            awk -v OFS="\\t" 'NR>1 && ${params.alt_col}+${params.ref_col}>0 {{print $1,$2+1,"-",$3,${params.alt_col},${params.ref_col}}}' | \
            gzip > {output}
        """


rule join_pileup_strands:
    input:
        expand(
            TEMPDIR / "pileup_by_strand/{{sample}}.{{reftype}}.{strand}.tsv.gz",
            strand=["fwd", "rev"],
        ),
    output:
        "report_sites/pileup/{sample}.{reftype}.tsv.gz",
    shell:
        """
        (echo -e "ref\tpos\tstrand\tmotif\tconverted\tunconverted" | gzip; cat {input} ) > {output}
        """


rule join_pileup_table:
    input:
        expand(
            "report_sites/pileup/{sample}.{{reftype}}.tsv.gz",
            sample=READS.keys(),
        ),
    output:
        "report_sites/joined/{reftype}.arrow",
    params:
        names=list(READS.keys()),
    threads: 8
    shell:
        """
        {PATH[joinPileup]} -f {input} -n {params.names} -o {output}
        """


# usage: filter_sites.py [-h] -f FILE -o OUTPUT [-u MIN_UNCON] [-d MIN_DEPTH] [-r MIN_RATIO] [-p MIN_PVAL]
rule filter_sites:
    input:
        "report_sites/joined/{reftype}.arrow",
    output:
        "report_sites/filtered/{reftype}.tsv.gz",
    params:
        min_uncon=config.get("cutoff", {}).get("min_uncon", 1),
        min_depth=config.get("cutoff", {}).get("min_depth", 10),
        min_ratio=config.get("cutoff", {}).get("min_ratio", 0.05),
        min_pval=config.get("cutoff", {}).get("min_pval", 1),
    threads: 8
    shell:
        """
        {PATH[filterSites]} -f {input} -o {output} -u {params.min_uncon} -d {params.min_depth} -r {params.min_ratio} -p {params.min_pval}
        """
