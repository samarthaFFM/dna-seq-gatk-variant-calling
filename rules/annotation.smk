rule snpeff:
    input:
        "filtered/all.vcf.gz",
    output:
        vcf=report("annotated/all.vcf.gz", caption="../report/vcf.rst", category="Calls"),
        csvstats="snpeff/all.csv"
    log:
        "logs/snpeff.log"
    params:
        reference=get_snpeff_database(),
        extra="-Xmx6g"
    resources:
        mem_gb=6
    wrapper:
        "0.27.1/bio/snpeff"


rule bcftools_view_vlr_for_annotation:
    input:
        "varlociraptor/filter/{filter}/{sample}-{unit}.fdr_{fdr}.{events}.{var_type}.bcf"
    output:
        "varlociraptor/filter/{filter}/{sample}-{unit}.fdr_{fdr}.{events}.{var_type}.vcf"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools view -o {output} -O v {input}"


rule snpeff_filtered_vlr:
    input:
        get_restrict_regions(),
        "varlociraptor/filter/{filter}/{sample}-{unit}.fdr_{fdr}.{events}.{var_type}.vcf"
    output:
        vcf=report("snpeff/varlociraptor/filter/{filter}/{sample}-{unit}.fdr_{fdr}.{events}.{var_type}.vcf.gz", caption="../report/vcf.rst", category="Calls"),
        csvstats="snpeff/varlociraptor/filter/{filter}/{sample}-{unit}.fdr_{fdr}.{events}.{var_type}.csv"
    log:
        "logs/snpeff/varlociraptor/filter/{filter}/{sample}-{unit}.fdr_{fdr}.{events}.{var_type}.snpeff.log"
    params:
        reference=get_snpeff_database(),
        extra="-Xmx6g -no-downstream -no-intergenic -no-intron -no-upstream -no-utr"
    resources:
        mem_gb=6
    wrapper:
        "0.27.1/bio/snpeff"


rule jannovar_filtered_vlr:
    input:
        vcf="varlociraptor/filter/{filter}/{sample}-{unit}.fdr_{fdr}.{events}.{var_type}.vcf"
    output:
        "jannovar/varlociraptor/filter/{filter}/{sample}-{unit}.fdr_{fdr}.{events}.{var_type}.vcf.gz"
    log:
        "logs/jannovar/varlociraptor/filter/{filter}/{sample}-{unit}.fdr_{fdr}.{events}.{var_type}.log"
    params:
        database=get_jannovar_database(), # path to jannovar reference dataset
        extra="--show-all"         # optional parameters
    wrapper:
        "0.35.0/bio/jannovar


rule igv_reports_vlr:
    input:
        vcf="{annotation}/varlociraptor/filter/{filter}/{sample}-{unit}.fdr_{fdr}.{events}.{var_type}.vcf.gz",
        bam="realigned/{sample}-{unit}.bam",
        ref=get_ref(),
        ref_idx=get_ref_idx(),
        regions=get_restrict_regions(),
        dbsnp=get_dbsnp()
    output:
        "{annotation}/varlociraptor/filter/{filter}/{sample}-{unit}.fdr_{fdr}.{events}.{var_type}.html"
    conda:
        "../envs/igv_reports.yaml"
    log:
        "logs/{annotation}/varlociraptor/filter/{filter}/{sample}-{unit}.fdr_{fdr}.{events}.{var_type}.igv_reports.log",
    shell:
        "( create_report {input.vcf}"
        "   {input.ref}"
        "   --flanking 1000"
        "   --sample-columns DP AF OBS"
        "   --info-columns PROB_ALT ANN"
        "   --tracks {input.bam} {input.regions} {input.dbsnp}"
        "   --output {output} ) 2>{log}"
        
