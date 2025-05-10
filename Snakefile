from collections import defaultdict
from snakemake.utils import validate


validate(workflow.configfiles[-1], "config.schema.yaml")


WORKDIR = os.path.relpath(
    config.get("workdir", "workspace"), os.path.dirname(workflow.configfiles[-1])
)


workdir: WORKDIR


def resolve_path(path):
    if os.path.isabs(path):
        return path
    else:
        path = os.path.expanduser(path)
        if os.path.isabs(path):
            return path
        else:
            return os.path.relpath(path, WORKDIR)


TEMPDIR = Path(
    os.path.relpath(
        config.get("tempdir", os.path.join(workflow.basedir, ".tmp")), workflow.basedir
    )
)
INTERNALDIR = Path("internal_files")

PATH = config["path"]

LIBTYPE = config["libtype"]
WITH_UMI = LIBTYPE in ["ECLIP10", "ECLIP6", "TAKARAV3", "SACSEQV3"]
MARKDUP = config.get("markdup", True)
# force markdup to be true for UMI
# if WITH_UMI:
#     MARKDUP = True
STRANDNESS = config.get("strandness", True)
GENE_NORC = config.get("gene_norc", True)
BASE_CHANGE = config.get("base_change", "A,G")
SPLICE_GENOME = config.get("splice_genome", True)
SPLICE_CONTAM = config.get("splice_contamination", False)


REF = config["reference"]
for k, v in REF.items():
    REF[k] = [resolve_path(v2) for v2 in v]

CONTAMINATION_FASTA = REF.get("contamination", [])
GENES_FASTA = REF.get("genes", [])
GENOME_FASTA = REF.get("genome", [])

REF_TYPES = (["contamination"] if len(CONTAMINATION_FASTA) > 0 else []) + [
    "genes",
    "genome",
]


SAMPLE2DATA = defaultdict(lambda: defaultdict(dict))
for s, v in config[f"samples"].items():
    for i, v2 in enumerate(v["data"], 1):
        r = f"run{i}"
        SAMPLE2DATA[str(s)][r] = {k: resolve_path(v3) for k, v3 in dict(v2).items()}


rule all:
    input:
        expand("report_reads/unmap/{sample}.count", sample=SAMPLE2DATA.keys()),
        expand(
            "report_reads/combined/{sample}.{reftype}.count",
            sample=SAMPLE2DATA.keys(),
            reftype=REF_TYPES,
        ),
        expand(
            "report_reads/dedup/{sample}.{reftype}.count",
            sample=SAMPLE2DATA.keys(),
            reftype=REF_TYPES,
        ),


rule cutadapt_SE:
    input:
        lambda wildcards: SAMPLE2DATA[wildcards.sample][wildcards.rn].get("R1", "/"),
    output:
        c=temp(TEMPDIR / "trimmed_reads/SE/{sample}_{rn}_R1.fq.gz"),
        s=INTERNALDIR / "discarded_reads/{sample}_{rn}_tooshort_R1.fq.gz",
        report="report_reads/trimming/{sample}_{rn}.json",
    params:
        library=LIBTYPE,
    threads: 16
    shell:
        """
        {PATH[cutseq]} -t {threads} -A {params.library:q} -m 20 --auto-rc -o {output.c} -s {output.s} --json-file {output.report} {input} 
        """


rule cutadapt_PE:
    input:
        lambda wildcards: SAMPLE2DATA[wildcards.sample][wildcards.rn].get("R1", "/"),
        lambda wildcards: SAMPLE2DATA[wildcards.sample][wildcards.rn].get("R2", "/"),
    output:
        c=[
            temp(TEMPDIR / "trimmed_reads/PE/{sample}_{rn}_R1.fq.gz"),
            temp(TEMPDIR / "trimmed_reads/PE/{sample}_{rn}_R2.fq.gz"),
        ],
        s=[
            INTERNALDIR / "discarded_reads/{sample}_{rn}_tooshort_R1.fq.gz",
            INTERNALDIR / "discarded_reads/{sample}_{rn}_tooshort_R2.fq.gz",
        ],
        report="report_reads/trimming/{sample}_{rn}.json",
    params:
        library=LIBTYPE,
    threads: 16
    shell:
        """
        {PATH[cutseq]} -t {threads} -A {params.library:q} -m 20 --auto-rc -o {output.c} -s {output.s} --json-file {output.report} {input} 
        """


# map to contamination (drop reads)


rule combine_contamination_fa:
    input:
        CONTAMINATION_FASTA,
    output:
        fa=INTERNALDIR / "reference_file/contamination.fa",
        fai=INTERNALDIR / "reference_file/contamination.fa.fai",
    shell:
        """
        zcat -f {input} > {output.fa}
        {PATH[samtools]} faidx {output.fa} --fai-idx {output.fai}
        """


rule build_contamination_hisat3n_index:
    input:
        INTERNALDIR / "reference_file/contamination.fa",
    output:
        INTERNALDIR / "reference_file/contamination.3n.GA.1.ht2",
    params:
        prefix=str(INTERNALDIR / "reference_file/contamination"),
    threads: 4
    shell:
        """
        rm -f {params.prefix}*.ht2
        {PATH[hisat3nbuild]} -p {threads} --base-change {BASE_CHANGE} {input} {params.prefix}
        """


rule hisat2_3n_mapping_contamination_SE:
    input:
        fq=TEMPDIR / "trimmed_reads/SE/{sample}_{rn}_R1.fq.gz",
        idx=INTERNALDIR / "reference_file/contamination.3n.GA.1.ht2",
    output:
        mapped=temp(
            TEMPDIR / "contamination_mapping/SE/{sample}_{rn}.contamination.bam"
        ),
        unmapped=temp(TEMPDIR / "contamination_mapping/SE/{sample}_{rn}.unmap.bam"),
        summary="report_reads/premap/{sample}_{rn}_SE.summary",
    params:
        index=str(INTERNALDIR / "reference_file/contamination"),
        directional=lambda wildcards: (
            "" if LIBTYPE == "UNSTRANDED" else "--directional-mapping"
        ),
        splice_args=(
            "--pen-noncansplice 20 --min-intronlen 20 --max-intronlen 20"
            if SPLICE_CONTAM
            else "--no-spliced-alignment"
        ),
    threads: 24
    shell:
        """
        {PATH[hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -U {input.fq} --base-change {BASE_CHANGE} {params.directional} {params.splice_args} \
            --np 0 --rdg 5,3 --rfg 5,3 --sp 9,3 --mp 3,1 --score-min L,-2,-0.8 |\
            {PATH[samtools]} view -@ {threads} -e '!flag.unmap && qlen-sclen >= 30 && [XM] * 15 < qlen-sclen' -O BAM -U {output.unmapped} -o {output.mapped}
        """


rule hisat2_3n_mapping_contamination_PE:
    input:
        fq1=TEMPDIR / "trimmed_reads/PE/{sample}_{rn}_R1.fq.gz",
        fq2=TEMPDIR / "trimmed_reads/PE/{sample}_{rn}_R2.fq.gz",
        idx=INTERNALDIR / "reference_file/contamination.3n.GA.1.ht2",
    output:
        mapped=temp(
            TEMPDIR / "contamination_mapping/PE/{sample}_{rn}.contamination.bam"
        ),
        unmapped=temp(TEMPDIR / "contamination_mapping/PE/{sample}_{rn}.unmap.bam"),
        summary="report_reads/premap/{sample}_{rn}_PE.summary",
    params:
        index=str(INTERNALDIR / "reference_file/contamination"),
        directional=lambda wildcards: (
            "" if LIBTYPE == "UNSTRANDED" else "--directional-mapping"
        ),
        splice_args=(
            "--pen-noncansplice 20 --min-intronlen 20 --max-intronlen 20"
            if SPLICE_CONTAM
            else "--no-spliced-alignment"
        ),
    threads: 24
    shell:
        """
        {PATH[hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -1 {input.fq1} -2 {input.fq2} --base-change {BASE_CHANGE} {params.directional} {params.splice_args} \
            --no-discordant --no-mixed \
            --np 0 --rdg 5,3 --rfg 5,3 --sp 9,3 --mp 3,1 --score-min L,-2,-0.8 |\
            {PATH[samtools]} view -@ {threads} -e 'flag.proper_pair && !flag.unmap && !flag.munmap && qlen-sclen >= 30 && [XM] * 15 < (qlen-sclen)' -O BAM -U {output.unmapped} -o {output.mapped}
        """


rule fix_hisat_tlen_bug_premap_SE:
    input:
        TEMPDIR / "contamination_mapping/SE/{sample}_{rn}.contamination.bam",
    output:
        temp(TEMPDIR / "contamination_fixmate/SE/{sample}_{rn}.bam"),
    threads: 4
    shell:
        "cp {input} {output}"


rule fix_hisat_tlen_bug_premap_PE:
    input:
        TEMPDIR / "contamination_mapping/PE/{sample}_{rn}.contamination.bam",
    output:
        temp(TEMPDIR / "contamination_fixmate/PE/{sample}_{rn}.bam"),
    threads: 4
    shell:
        "{PATH[samtools]} fixmate -@ {threads} -m -O BAM {input} {output}"


rule sort_bam_contamination:
    input:
        TEMPDIR / "contamination_fixmate/{mode}/{sample}_{rn}.bam",
    output:
        temp(TEMPDIR / "contamination_bam/{mode}/{sample}_{rn}.bam"),
    threads: 24
    shell:
        """
        {PATH[samtools]} sort -@ {threads} -m 3G -O BAM -o {output} {input}
        """


rule extract_unmapped_from_premap_SE:
    input:
        un=TEMPDIR / "contamination_mapping/SE/{sample}_{rn}.unmap.bam",
    output:
        r1=temp(TEMPDIR / "prefilter_reads/SE/{sample}_{rn}_R1.fq.gz"),
    shell:
        """
        {PATH[samtools]} fastq -0 {output.r1} -n {input}
        """


rule extract_unmapped_from_premap_PE:
    input:
        un=TEMPDIR / "contamination_mapping/PE/{sample}_{rn}.unmap.bam",
    output:
        r1=temp(TEMPDIR / "prefilter_reads/PE/{sample}_{rn}_R1.fq.gz"),
        r2=temp(TEMPDIR / "prefilter_reads/PE/{sample}_{rn}_R2.fq.gz"),
    shell:
        """
        {PATH[samtools]} fastq -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n {input}
        """


# map to genes


rule prepare_genes_index:
    input:
        GENES_FASTA,
    output:
        fa=INTERNALDIR / "genes_index/genes.fa",
        index=INTERNALDIR / "genes_index/genes.3n.CT.1.ht2",
    params:
        prefix=INTERNALDIR / "genes_index/genes",
    threads: 4
    shell:
        """
        zcat -f {input} | sed 's/(//g' | sed 's/)//g' > {output.fa}
        {PATH[samtools]} faidx {output.fa}
        rm -f {params.prefix}.3n.*.ht2
        {PATH[hisat3nbuild]} -p {threads} --base-change {BASE_CHANGE} {output.fa} {params.prefix}
        """


rule hisat2_3n_mapping_genes_SE:
    input:
        fq1=(
            TEMPDIR / "prefilter_reads/SE/{sample}_{rn}_R1.fq.gz"
            if len(CONTAMINATION_FASTA) > 0
            else TEMPDIR / "trimmed_reads/SE/{sample}_{rn}_R1.fq.gz"
        ),
        index=INTERNALDIR / "genes_index/genes.3n.CT.1.ht2",
    output:
        bam=temp(TEMPDIR / "genes_mapping/SE/{sample}_{rn}.bam"),
        summary="report_reads/mapping/{sample}_{rn}.genes.summary",
    params:
        index=INTERNALDIR / "genes_index/genes",
        directional="--directional-mapping" if STRANDNESS else "",
    threads: 24
    shell:
        """
        {PATH[hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -U {input.fq1} --base-change {BASE_CHANGE} {params.directional} \
            --avoid-pseudogene --no-softclip --no-spliced-alignment --np 0 --rdg 5,3 --rfg 5,3 --mp 3,1 --score-min L,-3,-0.75 |\
            {PATH[samtools]} view -@ {threads} -O BAM -o {output.bam}
        """


rule hisat2_3n_mapping_genes_PE:
    input:
        fq1=(
            TEMPDIR / "prefilter_reads/PE/{sample}_{rn}_R1.fq.gz"
            if len(CONTAMINATION_FASTA) > 0
            else TEMPDIR / "trimmed_reads/PE/{sample}_{rn}_R1.fq.gz"
        ),
        fq2=(
            TEMPDIR / "prefilter_reads/PE/{sample}_{rn}_R2.fq.gz"
            if len(CONTAMINATION_FASTA) > 0
            else TEMPDIR / "trimmed_reads/PE/{sample}_{rn}_R2.fq.gz"
        ),
        index=INTERNALDIR / "genes_index/genes.3n.CT.1.ht2",
    output:
        bam=temp(TEMPDIR / "genes_mapping/PE/{sample}_{rn}.bam"),
        summary="report_reads/mapping/{sample}_{rn}.genes.summary",
    params:
        index=INTERNALDIR / "genes_index/genes",
        directional="--directional-mapping" if STRANDNESS else "",
    threads: 24
    shell:
        """
        {PATH[hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -1 {input.fq1} -2 {input.fq2} --base-change {BASE_CHANGE} {params.directional} \
            --avoid-pseudogene --no-softclip --no-spliced-alignment --np 0 --rdg 5,3 --rfg 5,3 --mp 3,1 --score-min L,-3,-0.75 |\
            {PATH[samtools]} view -@ {threads} -O BAM -o {output.bam}
        """


rule filter_and_sort_bam_genes:
    input:
        TEMPDIR / "genes_mapping/{mode}/{sample}_{rn}.bam",
    output:
        unmap=temp(TEMPDIR / "genes_unmapped/{mode}/{sample}_{rn}.bam"),
        mapped=temp(TEMPDIR / "genes_bam/{mode}/{sample}_{rn}.bam"),
    params:
        flt=lambda wildcards: (
            "flag.proper_pair && !flag.unmap && !flag.munmap"
            + (
                " && !flag.read1 != !flag.reverse"
                if (GENE_NORC and STRANDNESS)
                else ""
            )
            if wildcards.mode == "PE"
            else "!flag.unmap"
            + (" && !flag.reverse" if (GENE_NORC and STRANDNESS) else "")
        ),
    threads: 16
    shell:
        """
        {PATH[samtools]} view -e '{params.flt}' -@ {threads} -u -U {output.unmap} {input} |\
            {PATH[samtools]} fixmate -@ {threads} -m -u - - |\
            {PATH[samtools]} sort -m 3G -O BAM -o {output.mapped}
        """


rule extract_unmapped_reads_SE_from_genes:
    input:
        TEMPDIR / "genes_unmapped/SE/{sample}_{rn}.bam",
    output:
        temp(TEMPDIR / "genes_unmapped/SE/{sample}_{rn}_R1.fq.gz"),
    threads: 4
    shell:
        """
        {PATH[samtools]} fastq -@ {threads} -0 {output} {input}
        """


rule extract_unmapped_reads_PE_from_genes:
    input:
        TEMPDIR / "genes_unmapped/PE/{sample}_{rn}.bam",
    output:
        r1=temp(TEMPDIR / "genes_unmapped/PE/{sample}_{rn}_R1.fq.gz"),
        r2=temp(TEMPDIR / "genes_unmapped/PE/{sample}_{rn}_R2.fq.gz"),
    threads: 4
    shell:
        """
        {PATH[samtools]} fastq -@ {threads} -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n {input}
        """


# map to genome


rule hisat2_3n_mapping_genome_SE:
    input:
        fq1=TEMPDIR / "genes_unmapped/SE/{sample}_{rn}_R1.fq.gz",
    output:
        bam=temp(TEMPDIR / "genome_mapping/SE/{sample}_{rn}.bam"),
        summary="report_reads/mapping/{sample}_{rn}_genome.summary",
    params:
        index=config.get("genome_index"),
        directional="--directional-mapping" if STRANDNESS else "",
        splice_args=(
            "--pen-noncansplice 20 --min-intronlen 20 --max-intronlen 20"
            if SPLICE_GENOME
            else "--no-spliced-alignment"
        ),
    threads: 24
    shell:
        """
        {PATH[hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -U {input.fq1} --base-change {BASE_CHANGE} {params.directional} {params.splice_args} \
            --avoid-pseudogene --np 0 --rdg 5,3 --rfg 5,3 --sp 9,3 --mp 3,1 --score-min L,-3,-0.75 |\
            {PATH[samtools]} view -@ {threads} -O BAM -o {output.bam}
        """


rule hisat2_3n_mapping_genome_PE:
    input:
        fq1=TEMPDIR / "genes_unmapped/PE/{sample}_{rn}_R1.fq.gz",
        fq2=TEMPDIR / "genes_unmapped/PE/{sample}_{rn}_R2.fq.gz",
    output:
        bam=temp(TEMPDIR / "genome_mapping/PE/{sample}_{rn}.bam"),
        summary="report_reads/mapping/{sample}_{rn}_genome.summary",
    params:
        index=config.get("genome_index"),
        directional="--directional-mapping" if STRANDNESS else "",
        splice_args=(
            "--pen-noncansplice 20 --min-intronlen 20 --max-intronlen 20"
            if SPLICE_GENOME
            else "--no-spliced-alignment"
        ),
    threads: 24
    shell:
        """
        {PATH[hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -1 {input.fq1} -2 {input.fq2} --base-change {BASE_CHANGE} {params.directional} {params.splice_args} \
            --avoid-pseudogene --np 0 --rdg 5,3 --rfg 5,3 --sp 9,3 --mp 3,1 --score-min L,-3,-0.75 |\
            {PATH[samtools]} view -@ {threads} -O BAM -o {output.bam}
        """


rule filter_and_sort_bam_genome:
    input:
        TEMPDIR / "genome_mapping/{mode}/{sample}_{rn}.bam",
    output:
        unmap=temp(TEMPDIR / "genome_unmapped/{mode}/{sample}_{rn}.bam"),
        mapped=temp(TEMPDIR / "genome_bam/{mode}/{sample}_{rn}.bam"),
    params:
        flt=lambda wildcards: (
            "flag.proper_pair && !flag.unmap && !flag.munmap"
            if wildcards.mode == "PE"
            else "!flag.unmap"
        ),
    threads: 16
    shell:
        # samtools collate -@ 4 -O -u example.bam | 
        """
        {PATH[samtools]} view -e '{params.flt}' -@ {threads} -u -U {output.unmap} {input} |\
            {PATH[samtools]} fixmate -@ {threads} -m -u - - |\
            {PATH[samtools]} sort -m 3G -O BAM -o {output.mapped}
        """


ruleorder: extract_unmapped_reads_PE_from_genome > extract_unmapped_reads_SE_from_genome


rule extract_unmapped_reads_SE_from_genome:
    input:
        TEMPDIR / "genome_unmapped/SE/{sample}_{rn}.bam",
    output:
        INTERNALDIR / "discarded_reads/{sample}_{rn}_unmap_R1.fq.gz",
    threads: 4
    shell:
        """
        {PATH[samtools]} fastq -@ {threads} -0 {output} {input}
        """


rule extract_unmapped_reads_PE_from_genome:
    input:
        TEMPDIR / "genome_unmapped/PE/{sample}_{rn}.bam",
    output:
        r1=INTERNALDIR / "discarded_reads/{sample}_{rn}_unmap_R1.fq.gz",
        r2=INTERNALDIR / "discarded_reads/{sample}_{rn}_unmap_R2.fq.gz",
    threads: 4
    shell:
        """
        {PATH[samtools]} fastq -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n {input}
        """


rule stat_unmapped:
    input:
        lambda wildcards: [
            INTERNALDIR / f"discarded_reads/{wildcards.sample}_{rn}_unmap_R1.fq.gz"
            for rn in SAMPLE2DATA[wildcards.sample].keys()
        ],
    output:
        "report_reads/unmap/{sample}.count",
    threads: 2
    shell:
        """
        zcat {input} | {PATH[samtools]} view -@ {threads} -c > {output}
        """


rule combine_runs:
    input:
        lambda wildcards: [
            (
                TEMPDIR / f"{wildcards.reftype}_bam/PE/{wildcards.sample}_{rn}.bam"
                if len(rs) == 2
                else TEMPDIR
                / f"{wildcards.reftype}_bam/SE/{wildcards.sample}_{rn}.bam"
            )
            for rn, rs in SAMPLE2DATA[wildcards.sample].items()
        ],
    output:
        temp(TEMPDIR / "combined_runs/{sample}.{reftype}.bam"),
    threads: 16
    shell:
        """
        {PATH[samtools]} merge -@ {threads} -f --write-index -o {output}##idx##{output}.bai {input}
        """


rule stat_combined:
    input:
        bam=TEMPDIR / "combined_runs/{sample}.{reftype}.bam",
    output:
        stat="report_reads/combined/{sample}.{reftype}.txt",
        n="report_reads/combined/{sample}.{reftype}.count",
    threads: 2
    shell:
        """
        {PATH[samtools]} flagstat -@ {threads} -O TSV {input} > {output.stat}
        {PATH[samtools]} view -@ {threads} -c -F 384 {input} > {output.n}
        """


rule drop_duplicates:
    input:
        bam=TEMPDIR / "combined_runs/{sample}.{reftype}.bam",
    output:
        bam=INTERNALDIR / "aligned_bam/{sample}.{reftype}.bam",
        txt="report_reads/dedup/{sample}.{reftype}.log",
    params:
        umi="--barcode-rgx '.*_([!-?A-~]+)'",
    threads: 16
    shell:
        # use bash, if MARKDUP is true, use samtools markdup, otherwise use cp
        """
        if [[ "{MARKDUP}" == "True" ]]; then
            {PATH[samtools]} markdup -@ {threads} -r -S --mode t --include-fails --duplicate-count --write-index {params.umi} -f {output.txt} {input} {output.bam}
        else
            cp {input.bam} {output.bam} && touch {output.txt}
        fi
        """


rule stat_dedup:
    input:
        INTERNALDIR / "aligned_bam/{sample}.{reftype}.bam",
    output:
        stat="report_reads/dedup/{sample}.{reftype}.txt",
        n="report_reads/dedup/{sample}.{reftype}.count",
    threads: 2
    shell:
        """
        {PATH[samtools]} flagstat -@ {threads} -O TSV {input} > {output.stat}
        {PATH[samtools]} view -@ {threads} -c -F 384 {input} > {output.n}
        """
