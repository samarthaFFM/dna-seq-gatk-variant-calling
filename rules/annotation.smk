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

rule download_jannovar_database:
    output:
        directory("data/ref/jannovar/{jannovar_db}")
    log:
        "logs/jannovar/{jannovar_db}/download.log"
    conda:
        "../envs/jannovar.yaml"
    resources:
        mem_gb=12
    shell:
        "jannovar download  -Xmx{resources.mem_gb}g -d {wildcards.jannovar_db} --download-dir {output}"
        

rule jannovar_filtered_vlr:
    input:
        vcf="varlociraptor/filter/{filter}/{sample}-{unit}.fdr_{fdr}.{events}.{var_type}.vcf",
        database="data/ref/jannovar/" + get_jannovar_database()
    output:
        "jannovar/varlociraptor/filter/{filter}/{sample}-{unit}.fdr_{fdr}.{events}.{var_type}.vcf.gz"
    log:
        "logs/jannovar/varlociraptor/filter/{filter}/{sample}-{unit}.fdr_{fdr}.{events}.{var_type}.log"
    params:
        database= lambda w, input: f"{input.database}/{get_jannovar_database().replace('/','_')}.ser", # path to jannovar reference dataset
        extra= lambda w, output: f"--show-all -Xms1g -Xmx6g"       # optional parameters
    resources:
        mem_gb=6
    wrapper:
        "0.35.0/bio/jannovar"


rule jannovar_filtered_vlr_fix_ANN_header:
    input:
        "jannovar/varlociraptor/filter/{filter}/{sample}-{unit}.fdr_{fdr}.{events}.{var_type}.vcf.gz"
    output:
        "jannovar/varlociraptor/filter/{filter}/{sample}-{unit}.fdr_{fdr}.{events}.{var_type}.fixed_ANN_header.vcf.gz"
    log:
        "logs/jannovar/varlociraptor/filter/{filter}/{sample}-{unit}.fdr_{fdr}.{events}.{var_type}.fixed_ANN_header.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "( bcftools view {input} |"
        "  sed -e 's/^##INFO=<ID=ANN,Number=1,Type/##INFO=<ID=ANN,Number=.,Type/' -e 's/^##FORMAT=<ID=DP,Number=1,Type/##FORMAT=<ID=DP,Number=A,Type/' |"
        "  bcftools view -O z -o {output} ) 2>{log}"


def parse_events(wildcards):
    events=wildcards.events.split('-')
    events.insert(0, "")
    return " PROB_".join( events )

rule igv_reports_vlr:
    input:
        vcf="{annotation}/varlociraptor/filter/{filter}/{sample}-{unit}.fdr_{fdr}.{events}.{var_type,[^\.]+}{fixed,(|\.fixed_ANN_header)}.vcf.gz",
        bam="realigned/{sample}-{unit}.bam",
        ref=get_ref(),
        ref_idx=get_ref_idx(),
        regions=get_restrict_regions(),
        dbsnp=get_dbsnp()
    output:
        "{annotation}/varlociraptor/filter/{filter}/{sample}-{unit}.fdr_{fdr}.{events}.{var_type}{fixed}.html"
    params:
        events=parse_events
    conda:
        "../envs/igv_reports.yaml"
    log:
        "logs/{annotation}/varlociraptor/filter/{filter}/{sample}-{unit}.fdr_{fdr}.{events}.{var_type}{fixed}.igv_reports.log",
    shell:
        "( create_report {input.vcf}"
        "   {input.ref}"
        "   --flanking 1000"
        "   --sample-columns DP AF OBS"
        "   --info-columns {params.events} ANN"
        "   --tracks {input.bam} {input.regions} {input.dbsnp}"
        "   --output {output} ) 2>{log}"
        
