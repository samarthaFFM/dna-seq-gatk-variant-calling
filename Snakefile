include: "rules/common.smk"

##### Target rules #####

rule all:
    input:
        "annotated/all.vcf.gz",
        "qc/multiqc.html",
        "qc/fastqc.html",
        "plots/depths.svg",
        "plots/allele-freqs.svg"


##### Modules #####

include: "rules/download.smk"
include: "rules/mapping.smk"
include: "rules/calling.smk"
include: "rules/filtering.smk"
include: "rules/stats.smk"
include: "rules/qc.smk"
include: "rules/annotation.smk"
