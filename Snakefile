from collections import defaultdict
from snakemake.utils import validate
from workflow_utils import preprocess_config


configfile: "default.yaml"


validate(workflow.configfiles[0], "config.schema.yaml")


TEMPDIR = Path(
    os.path.relpath(
        config.get("tempdir", os.path.join(workflow.basedir, ".tmp")), workflow.basedir
    )
)
INTERNALDIR = Path("internal_files")


if workflow._workdir_handler is None:

    WORKDIR = os.path.relpath(
        config.get("workdir", "workspace"), os.path.dirname(workflow.configfiles[0])
    )

    workdir: WORKDIR

else:
    WORKDIR = workflow._workdir_handler.workdir


preprocess_config(config, WORKDIR)


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
REF = config.get("_REF", {})
READS = config.get("_READS", {})


rule all:
    input:
        expand("report_reads/unmap/{sample}.count", sample=READS.keys()),
        expand(
            "report_reads/combined/{sample}.{reftype}.count",
            sample=READS.keys(),
            reftype=REF.keys(),
        ),
        expand(
            "report_reads/dedup/{sample}.{reftype}.count",
            sample=READS.keys(),
            reftype=REF.keys(),
        ),


rule cutadapt_SE:
    input:
        lambda wildcards: READS[wildcards.sample][wildcards.rn].get("R1", "/"),
    output:
        c=temp(TEMPDIR / "trimmed_reads/SE/{sample}_{rn}_R1.fq.gz"),
        s=INTERNALDIR / "discarded_reads/{sample}_{rn}_tooshort_R1.fq.gz",
        report="report_reads/trimming/{sample}_{rn}.json",
    params:
        library=LIBTYPE,
    threads: 16
    shell:
        """
        {config[path][cutseq]} -t {threads} -A {params.library:q} -m 20 --auto-rc -o {output.c} -s {output.s} --json-file {output.report} {input} 
        """


rule cutadapt_PE:
    input:
        lambda wildcards: READS[wildcards.sample][wildcards.rn].get("R1", "/"),
        lambda wildcards: READS[wildcards.sample][wildcards.rn].get("R2", "/"),
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
        {config[path][cutseq]} -t {threads} -A {params.library:q} -m 20 --auto-rc -o {output.c} -s {output.s} --json-file {output.report} {input} 
        """


# map to contamination (drop reads)


rule combine_contamination_fa:
    input:
        REF.get("contamination", []),
    output:
        fa=INTERNALDIR / "reference_file/contamination.fa",
        fai=INTERNALDIR / "reference_file/contamination.fa.fai",
    shell:
        """
        zcat -f {input} > {output.fa}
        {config[path][samtools]} faidx {output.fa} --fai-idx {output.fai}
        """


rule build_contamination_hisat3n_index:
    input:
        INTERNALDIR / "reference_file/contamination.fa",
    output:
        INTERNALDIR / "reference_file/contamination.3n.GA.1.ht2",
    params:
        prefix=str(INTERNALDIR / "reference_file/contamination"),
        base_change=BASE_CHANGE,
    threads: 4
    shell:
        """
        rm -f {params.prefix}*.ht2
        {config[path][hisat3nbuild]} -p {threads} --base-change {params.base_change} {input} {params.prefix}
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
        base_change=BASE_CHANGE,
    threads: 24
    shell:
        """
        {config[path][hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -U {input.fq} --base-change {params.base_change} {params.directional} {params.splice_args} \
            --np 0 --rdg 5,3 --rfg 5,3 --sp 9,3 --mp 3,1 --score-min L,-2,-0.8 |\
            {config[path][samtools]} view -@ {threads} -e '!flag.unmap && qlen-sclen >= 30 && [XM] * 15 < qlen-sclen' -O BAM -U {output.unmapped} -o {output.mapped}
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
        base_change=BASE_CHANGE,
    threads: 24
    shell:
        """
        {config[path][hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -1 {input.fq1} -2 {input.fq2} --base-change {params.base_change} {params.directional} {params.splice_args} \
            --no-discordant --no-mixed \
            --np 0 --rdg 5,3 --rfg 5,3 --sp 9,3 --mp 3,1 --score-min L,-2,-0.8 |\
            {config[path][samtools]} view -@ {threads} -e 'flag.proper_pair && !flag.unmap && !flag.munmap && qlen-sclen >= 30 && [XM] * 15 < (qlen-sclen)' -O BAM -U {output.unmapped} -o {output.mapped}
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
        "{config[path][samtools]} fixmate -@ {threads} -m -O BAM {input} {output}"


rule sort_bam_contamination:
    input:
        TEMPDIR / "contamination_fixmate/{mode}/{sample}_{rn}.bam",
    output:
        temp(TEMPDIR / "contamination_bam/{mode}/{sample}_{rn}.bam"),
    threads: 24
    shell:
        """
        {config[path][samtools]} sort -@ {threads} -m 3G -O BAM -o {output} {input}
        """


rule extract_unmapped_from_premap_SE:
    input:
        un=TEMPDIR / "contamination_mapping/SE/{sample}_{rn}.unmap.bam",
    output:
        r1=temp(TEMPDIR / "prefilter_reads/SE/{sample}_{rn}_R1.fq.gz"),
    shell:
        """
        {config[path][samtools]} fastq -0 {output.r1} -n {input}
        """


rule extract_unmapped_from_premap_PE:
    input:
        un=TEMPDIR / "contamination_mapping/PE/{sample}_{rn}.unmap.bam",
    output:
        r1=temp(TEMPDIR / "prefilter_reads/PE/{sample}_{rn}_R1.fq.gz"),
        r2=temp(TEMPDIR / "prefilter_reads/PE/{sample}_{rn}_R2.fq.gz"),
    shell:
        """
        {config[path][samtools]} fastq -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n {input}
        """


# map to genes


rule prepare_genes_index:
    input:
        REF.get("genes", []),
    output:
        fa=INTERNALDIR / "genes_index/genes.fa",
        index=INTERNALDIR / "genes_index/genes.3n.CT.1.ht2",
    params:
        prefix=INTERNALDIR / "genes_index/genes",
        base_change=BASE_CHANGE,
    threads: 4
    shell:
        """
        zcat -f {input} | sed 's/(//g' | sed 's/)//g' > {output.fa}
        {config[path][samtools]} faidx {output.fa}
        rm -f {params.prefix}.3n.*.ht2
        {config[path][hisat3nbuild]} -p {threads} --base-change {params.base_change} {output.fa} {params.prefix}
        """


rule hisat2_3n_mapping_genes_SE:
    input:
        fq1=(
            TEMPDIR / "prefilter_reads/SE/{sample}_{rn}_R1.fq.gz"
            if "contamination" in REF
            else TEMPDIR / "trimmed_reads/SE/{sample}_{rn}_R1.fq.gz"
        ),
        index=INTERNALDIR / "genes_index/genes.3n.CT.1.ht2",
    output:
        bam=temp(TEMPDIR / "genes_mapping/SE/{sample}_{rn}.bam"),
        summary="report_reads/mapping/{sample}_{rn}.genes.summary",
    params:
        index=INTERNALDIR / "genes_index/genes",
        directional="--directional-mapping" if STRANDNESS else "",
        base_change=BASE_CHANGE,
    threads: 24
    shell:
        """
        {config[path][hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -U {input.fq1} --base-change {params.base_change} {params.directional} \
            --avoid-pseudogene --no-softclip --no-spliced-alignment --np 0 --rdg 5,3 --rfg 5,3 --mp 3,1 --score-min L,-3,-0.75 |\
            {config[path][samtools]} view -@ {threads} -O BAM -o {output.bam}
        """


rule hisat2_3n_mapping_genes_PE:
    input:
        fq1=(
            TEMPDIR / "prefilter_reads/PE/{sample}_{rn}_R1.fq.gz"
            if "contamination" in REF
            else TEMPDIR / "trimmed_reads/PE/{sample}_{rn}_R1.fq.gz"
        ),
        fq2=(
            TEMPDIR / "prefilter_reads/PE/{sample}_{rn}_R2.fq.gz"
            if "contamination" in REF
            else TEMPDIR / "trimmed_reads/PE/{sample}_{rn}_R2.fq.gz"
        ),
        index=INTERNALDIR / "genes_index/genes.3n.CT.1.ht2",
    output:
        bam=temp(TEMPDIR / "genes_mapping/PE/{sample}_{rn}.bam"),
        summary="report_reads/mapping/{sample}_{rn}.genes.summary",
    params:
        index=INTERNALDIR / "genes_index/genes",
        directional="--directional-mapping" if STRANDNESS else "",
        base_change=BASE_CHANGE,
    threads: 24
    shell:
        """
        {config[path][hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -1 {input.fq1} -2 {input.fq2} --base-change {params.base_change} {params.directional} \
            --avoid-pseudogene --no-softclip --no-spliced-alignment --np 0 --rdg 5,3 --rfg 5,3 --mp 3,1 --score-min L,-3,-0.75 |\
            {config[path][samtools]} view -@ {threads} -O BAM -o {output.bam}
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
        {config[path][samtools]} view -e '{params.flt}' -@ {threads} -u -U {output.unmap} {input} |\
            {config[path][samtools]} fixmate -@ {threads} -m -u - - |\
            {config[path][samtools]} sort -m 3G -O BAM -o {output.mapped}
        """


rule extract_unmapped_reads_SE_from_genes:
    input:
        TEMPDIR / "genes_unmapped/SE/{sample}_{rn}.bam",
    output:
        temp(TEMPDIR / "genes_unmapped/SE/{sample}_{rn}_R1.fq.gz"),
    threads: 4
    shell:
        """
        {config[path][samtools]} fastq -@ {threads} -0 {output} {input}
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
        {config[path][samtools]} fastq -@ {threads} -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n {input}
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
        base_change=BASE_CHANGE,
    threads: 24
    shell:
        """
        {config[path][hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -U {input.fq1} --base-change {params.base_change} {params.directional} {params.splice_args} \
            --avoid-pseudogene --np 0 --rdg 5,3 --rfg 5,3 --sp 9,3 --mp 3,1 --score-min L,-3,-0.75 |\
            {config[path][samtools]} view -@ {threads} -O BAM -o {output.bam}
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
        base_change=BASE_CHANGE,
    threads: 24
    shell:
        """
        {config[path][hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -1 {input.fq1} -2 {input.fq2} --base-change {params.base_change} {params.directional} {params.splice_args} \
            --avoid-pseudogene --np 0 --rdg 5,3 --rfg 5,3 --sp 9,3 --mp 3,1 --score-min L,-3,-0.75 |\
            {config[path][samtools]} view -@ {threads} -O BAM -o {output.bam}
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
        {config[path][samtools]} view -e '{params.flt}' -@ {threads} -u -U {output.unmap} {input} |\
            {config[path][samtools]} fixmate -@ {threads} -m -u - - |\
            {config[path][samtools]} sort -m 3G -O BAM -o {output.mapped}
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
        {config[path][samtools]} fastq -@ {threads} -0 {output} {input}
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
        {config[path][samtools]} fastq -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n {input}
        """


rule stat_unmapped:
    input:
        lambda wildcards: [
            INTERNALDIR / f"discarded_reads/{wildcards.sample}_{rn}_unmap_R1.fq.gz"
            for rn in READS[wildcards.sample].keys()
        ],
    output:
        "report_reads/unmap/{sample}.count",
    threads: 2
    shell:
        """
        zcat {input} | {config[path][samtools]} view -@ {threads} -c > {output}
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
            for rn, rs in READS[wildcards.sample].items()
        ],
    output:
        bam=temp(TEMPDIR / "combined_runs/{sample}.{reftype}.bam"),
        bai=temp(TEMPDIR / "combined_runs/{sample}.{reftype}.bam.bai"),
    threads: 16
    shell:
        """
        {config[path][samtools]} merge -@ {threads} -f --write-index -o {output.bam}##idx##{output.bai} {input}
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
        {config[path][samtools]} flagstat -@ {threads} -O TSV {input} > {output.stat}
        {config[path][samtools]} view -@ {threads} -c -F 384 {input} > {output.n}
        """


rule drop_duplicates:
    input:
        bam=TEMPDIR / "combined_runs/{sample}.{reftype}.bam",
    output:
        bam=INTERNALDIR / "aligned_bam/{sample}.{reftype}.bam",
        txt="report_reads/dedup/{sample}.{reftype}.log",
    params:
        umi="--barcode-rgx '.*_([!-?A-~]+)'",
        markdup=MARKDUP,
    threads: 16
    shell:
        """
        if [[ "{params.markdup}" == "True" ]]; then
            {config[path][samtools]} markdup -@ {threads} -r -S --mode t --include-fails --duplicate-count --write-index {params.umi} -f {output.txt} {input} {output.bam}
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
        {config[path][samtools]} flagstat -@ {threads} -O TSV {input} > {output.stat}
        {config[path][samtools]} view -@ {threads} -c -F 384 {input} > {output.n}
        """
