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
        genome_path + "{genome}" + get_ref_ext()
    output:
        genome_path + "{genome}.amb",
        genome_path + "{genome}.ann",
        genome_path + "{genome}.bwt",
        genome_path + "{genome}.pac",
        genome_path + "{genome}.sa"
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
        ref_bwt = get_ref_basename() + ".bwt"
    output:
        temp("mapped/{sample}-{unit}.sorted.bam")
    log:
        "logs/bwa_mem/{sample}-{unit}.log"
    params:
        index = get_ref_basename(),
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


rule create_dict:
    input:
        genome_path + "{genome}" + get_ref_ext()
    output:
        genome_path + "{genome}.dict"
    log:
        "logs/picard_create_sequence_dictionary/{genome}.log"
    conda:
        "../envs/picard.yaml"
    shell:
        "picard CreateSequenceDictionary "
        "  R={input}"
        "  O={output}"


rule recalibrate_base_qualities:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref=get_ref(),
        dict=get_ref_basename() + ".dict",
        known=get_dbsnp(),
        known_idx=get_dbsnp_idx()
    output:
        bam=protected("recal/{sample}-{unit}.bam")
    params:
        extra=get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"],
        java_opts="-Xmx60G -XX:ParallelGCThreads=8"
    resources:
        mem_gb=60
    threads: 8
    log:
        "logs/gatk/bqsr/{sample}-{unit}.log"
    wrapper:
        "0.34.0/bio/gatk/baserecalibrator"


rule realigner_target_creator:
    input:
        bam="recal/{sample}-{unit}.bam",
        bai="recal/{sample}-{unit}.bai",
        ref=get_ref(),
        known=get_dbsnp(),
        known_idx=get_dbsnp_idx()
    output:
        targets="realigned/{sample}-{unit}.intervals"
    threads: 4
    params:
        threads=7
    log:
        "logs/gatk/realigner_target_creator/{sample}-{unit}.log"
    shell:
        "( gatk3 -T RealignerTargetCreator"
        "   -R {input.ref}"
        "   -I {input.bam}"
        "   -known {input.known}"
        "   -nt {params.threads}"
        "   --disable_auto_index_creation_and_locking_when_reading_rods"
        "   -o {output.targets} ) 2> {log}"


rule indel_realigner:
    input:
        bam="recal/{sample}-{unit}.bam",
        bai="recal/{sample}-{unit}.bai",
        ref=get_ref(),
        known=get_dbsnp(),
        known_idx=get_dbsnp_idx(),
        targets="realigned/{sample}-{unit}.intervals"
    output:
        bam="realigned/{sample}-{unit}.bam"
    log:
        "logs/gatk/indel_realigner/{sample}-{unit}.log"
    shell:
        "( gatk3 -T IndelRealigner"
        "   -R {input.ref}"
        "   -I {input.bam}"
        "   -known {input.known}"
        "   -targetIntervals {input.targets}"
        "   -o {output.bam} ) 2>{log}"




rule samtools_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    wrapper:
        "0.27.1/bio/samtools/index"
