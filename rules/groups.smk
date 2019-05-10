rule split_into_groups:
    input:
        "annotated/all.vcf.gz"
    output:
        "group/{group_id}.bcf"
    log:
        "bcftools_view/split_by_group/{group_id}.log"
    params:
        bcftools_extra_params = lambda w: 
            f"-O b --trim-alt-alleles --min-ac 1:nref --samples {','.join( samples.loc[samples['group'] == w.group_id, 'sample' ] )}"
    wrapper:
        "0.34.0/bio/bcftools/view"
