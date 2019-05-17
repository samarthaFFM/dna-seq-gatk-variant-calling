import pandas as pd
import os.path as path
from snakemake.utils import validate

report: "../report/workflow.rst"

###### Config file and sample sheets #####
configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep = '\t').set_index("sample", drop=False)
if 'group' in samples:
    samples.group = samples.group.astype(str)
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_csv(config["units"], dtype=str, sep = '\t').set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="../schemas/units.schema.yaml")

references = pd.read_csv(config["references"], sep = '\t').set_index("reference_type", drop=False)
validate(references, schema="../schemas/references.schema.yaml")


##### Wildcard constraints #####
wildcard_constraints:
    vartype="snvs|indels",
    sample="|".join(samples.index),
    unit="|".join(units["unit"]),
    reference_type="|".join(['genome', 'dbsnp']),
    fasta_ext="|".join(['fa', 'fasta', 'fna']),
    idx_ext="|".join(['fai', 'csi', 'tbi', 'idx'])


##### Helper functions #####

def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}

genome_path = "data/ref/genome/"

def get_ref():
    """Get the file name of the genome fasta file in references.tsv."""
    return genome_path + references.loc['genome'].get("file")

def get_ref_basename():
    """Get the basename of the genome fasta file in references.tsv without its extension."""
    return genome_path + path.splitext( references.loc['genome'].get("file") )[0]

def get_ref_ext():
    """Get the extension of the genome fasta file in references.tsv."""
    return path.splitext( references.loc['genome'].get("file") )[1]

def get_ref_idx():
    """Get the file name of the genome fasta index file in references.tsv."""
    return genome_path + references.loc['genome'].get("index")

def get_snpeff_database():
    """Get the snpeff reference database name for the genome fasta file in references.tsv."""
    return references.loc['genome'].get("snpeff_database")

   
dbsnp_path = "data/ref/dbsnp/"

def get_dbsnp():
    """Get the file name of the genome fasta file in references.tsv."""
    return dbsnp_path + references.loc['dbsnp'].get("file")

def get_dbsnp_idx():
    """Get the file name of the genome fasta file in references.tsv."""
    return dbsnp_path + references.loc['dbsnp'].get("index")


def get_restrict_regions_basename():
    """Get the basename of the genome fasta file in references.tsv without its extension."""
    return path.splitext( config["processing"]["restrict-regions"] )[0]


def is_single_end(sample, unit):
    """Return True if sample-unit is single end."""
    return pd.isnull(units.loc[(sample, unit), "fq2"])


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample,
        platform=units.loc[(wildcards.sample, wildcards.unit), "platform"])


def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand("trimmed/{sample}-{unit}.{group}.fastq.gz",
                      group=[1, 2], **wildcards)
    # single end sample
    return "trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)


def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand("recal/{sample}-{unit}.bam",
                  sample=wildcards.sample,
                  unit=units.loc[wildcards.sample].unit)


def get_regions_param(regions=config["processing"].get("restrict-regions"), default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        padding = config["processing"].get("region-padding")
        if padding:
            params += "--interval-padding {}".format(padding)
        return params
    return default


def get_call_variants_params(wildcards, input):
    config_params = config['params']['gatk']['HaplotypeCaller']
    regions_params = get_regions_param(regions=input.regions, default=f"--intervals {wildcards.contig}")
    return f"{config_params} {regions_params}"

def get_recal_input(bai=False):
    # case 1: no duplicate removal
    f = "mapped/{sample}-{unit}.sorted.bam"
    if config["processing"]["remove-duplicates"]:
        # case 2: remove duplicates
        f = "dedup/{sample}-{unit}.bam"
    if bai:
        if config["processing"].get("restrict-regions"):
            # case 3: need an index because random access is required
            f += ".bai"
            return f
        else:
            # case 4: no index needed
            return []
    else:
        return f
