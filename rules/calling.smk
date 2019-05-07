if "restrict-regions" in config["processing"]:
    rule compose_regions:
        input:
            regions = config["processing"]["restrict-regions"]
        output:
            "called/{contig}.regions.bed"
        conda:
            "../envs/bedops.yaml"
        shell:
            "bedextract {wildcards.contig} {input.regions} > {output}"


rule call_variants:
    input:
        bam=get_sample_bams,
        ref=get_ref(),
        known=get_dbsnp(),
        known_idx=get_dbsnp_idx(),
        regions="called/{contig}.regions.bed" if config["processing"].get("restrict-regions") else []
    output:
        gvcf=protected("called/{sample}.{contig}.g.vcf.gz")
    log:
        "logs/gatk/haplotypecaller/{sample}.{contig}.log"
    params:
        extra=get_call_variants_params
    wrapper:
        "0.27.1/bio/gatk/haplotypecaller"


rule combine_calls:
    input:
        ref=get_ref(),
        gvcfs=expand("called/{sample}.{{contig}}.g.vcf.gz", sample=samples.index)
    output:
        gvcf="called/all.{contig}.g.vcf.gz"
    log:
        "logs/gatk/combinegvcfs.{contig}.log"
    wrapper:
        "0.27.1/bio/gatk/combinegvcfs"


rule genotype_variants:
    input:
        ref=get_ref(),
        gvcf="called/all.{contig}.g.vcf.gz"
    output:
        vcf=temp("genotyped/all.{contig}.vcf.gz")
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"]
    log:
        "logs/gatk/genotypegvcfs.{contig}.log"
    wrapper:
        "0.27.1/bio/gatk/genotypegvcfs"


def generate_genotype_interval_vcfs(wildcards):
    """
    Get the contigs list from the genome fasta index file in references.tsv and generate vcf file names for intervals to genotype.
    """
    return expand("genotyped/all.{contig}.vcf.gz",
                    contig = pd.read_csv(checkpoints.idx_download.get(
                                reference_type="genome",
                                reference_file=references.loc['genome'].get("file"),
                                idx_ext=path.splitext( references.loc['genome'].get("index") )[1].lstrip('.')
                                ).output.idx,
                            header=None, usecols=[0], squeeze=True, dtype=str, sep = '\t')
                    )

rule merge_variants:
    input:
        vcf=generate_genotype_interval_vcfs
    output:
        vcf="genotyped/all.vcf.gz"
    log:
        "logs/picard/merge-genotyped.log"
    wrapper:
        "0.27.1/bio/picard/mergevcfs"
