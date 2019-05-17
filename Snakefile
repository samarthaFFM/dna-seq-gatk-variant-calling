include: "rules/common.smk"

##### Target rules #####

rule all:
    input:
        "annotated/all.vcf.gz",
        expand("group/{group_id}.bcf",  group_id = samples['group'].unique() ),
        expand("varlociraptor/filter/fdr/{u.sample}-{u.unit}.fdr_{fdr}.{events}.{var_type}.bcf",
                u = units.itertuples(),
                fdr = ["0-05", "0-001"],
                events = ["HET-HOM_ALT"],
                var_type = ["SNV", "INS", "DEL"] ),
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
