rule trim_reads_se:
    input:
        unpack(get_fastq)
    output:
        temp("trimmed/{sample}-{unit}.fastq.gz")
    params:
        extra="",
        **config["params"]["trimmomatic"]["se"]
    log:
        "logs/trimmomatic/{sample}-{unit}.log"
    wrapper:
        "0.30.0/bio/trimmomatic/se"


rule trim_reads_pe:
    input:
        unpack(get_fastq)
    output:
        r1=temp("trimmed/{sample}-{unit}.1.fastq.gz"),
        r2=temp("trimmed/{sample}-{unit}.2.fastq.gz"),
        r1_unpaired=temp("trimmed/{sample}-{unit}.1.unpaired.fastq.gz"),
        r2_unpaired=temp("trimmed/{sample}-{unit}.2.unpaired.fastq.gz"),
        trimlog="trimmed/{sample}-{unit}.trimlog.txt"
    params:
        extra=lambda w, output: "-trimlog {}".format(output.trimlog),
        **config["params"]["trimmomatic"]["pe"]
    log:
        "logs/trimmomatic/{sample}-{unit}.log"
    wrapper:
        "0.30.0/bio/trimmomatic/pe"

rule bwa_index:
    input:
        "data/ref/genome/{genome}.fasta"
    output:
        "data/ref/genome/{genome}.amb",
        "data/ref/genome/{genome}.ann",
        "data/ref/genome/{genome}.bwt",
        "data/ref/genome/{genome}.pac",
        "data/ref/genome/{genome}.sa"
    log:
        "logs/bwa_index/{genome}.log"
    params:
        prefix="data/ref/genome/{genome}",
        algorithm="bwtsw"
    wrapper:
        "0.33.0/bio/bwa/index"

rule map_reads:
    input:
        reads = get_trimmed_reads,
        ref = "data/ref/genome/" + get_ref
    output:
        temp("mapped/{sample}-{unit}.sorted.bam")
    log:
        "logs/bwa_mem/{sample}-{unit}.log"
    params:
        index=input.ref,
        extra=get_read_group,
        sort="samtools",
        sort_order="coordinate"
    threads: 8
    wrapper:
        "0.27.1/bio/bwa/mem"


rule mark_duplicates:
    input:
        "mapped/{sample}-{unit}.sorted.bam"
    output:
        bam=temp("dedup/{sample}-{unit}.bam"),
        metrics="qc/dedup/{sample}-{unit}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}-{unit}.log"
    params:
        config["params"]["picard"]["MarkDuplicates"]
    wrapper:
        "0.26.1/bio/picard/markduplicates"


rule recalibrate_base_qualities:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref="data/ref/genome/" + get_ref,
        known="data/ref/dbsnp/" + get_dbsnp
    output:
        bam=protected("recal/{sample}-{unit}.bam")
    params:
        extra=get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"]
    log:
        "logs/gatk/bqsr/{sample}-{unit}.log"
    wrapper:
        "0.27.1/bio/gatk/baserecalibrator"


rule samtools_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    wrapper:
        "0.27.1/bio/samtools/index"
