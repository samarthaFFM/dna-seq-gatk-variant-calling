include: "rules/common.smk"

##### Target rules #####

rule all:
    input:
        "annotated/all.vcf.gz",
        expand("group/{group_id}.bcf",  group_id = samples['group'].unique() ),
        expand("{annotation}/varlociraptor/filter/fdr/{u.sample}-{u.unit}.fdr_{fdr}.{events}.{var_type}.html",
                annotation = ["snpeff"],
                u = units.itertuples(),
                fdr = ["0-05", "0-001"],
                events = ["ALT-VERY_RARE"],
                var_type = ["SNV", "INS", "DEL"] ),
        expand("{annotation}/varlociraptor/filter/fdr/{u.sample}-{u.unit}.fdr_{fdr}.{events}.{var_type}.fixed_ANN_header.html",
                annotation = ["jannovar"],
                u = units.itertuples(),
                fdr = ["0-05", "0-001"],
                events = ["ALT-VERY_RARE"],
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
