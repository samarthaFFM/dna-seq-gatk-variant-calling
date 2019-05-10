include: "rules/common.smk"

##### Target rules #####

rule all:
    input:
        "annotated/all.vcf.gz",
        expand("group/{group_id}.bcf",  = samples['group'].unique() ),
        "qc/multiqc.html",
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
include: "rules/groups.smk"
