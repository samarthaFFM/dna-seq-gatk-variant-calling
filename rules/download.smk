import os.path as path

#No download; keep direct links on sample file
# rule(s) to download data from ENA based on the samples.tsv file

# rule ena_download:
#     output:
#         "data/raw/{sample}.1.fq.gz",
#         "data/raw/{sample}.2.fq.gz"
#     conda:
#         "../envs/curl.yaml"
#    params:
#        fq1 = lambda wildcards: samples.loc[wildcards.sample].get("fastq_ftp_1"),
#         fq2 = lambda wildcards: samples.loc[wildcards.sample].get("fastq_ftp_2")
#     resources:
#         parallel_download_connections = 1
#     log:
#         "logs/ena_download/{sample}.log"
#     shell:
#         """
#         ( curl -sS --output {output[0]} {params.fq1} ) 2> {log};
#         ( curl -sS --output {output[1]} {params.fq2} ) 2>> {log};
#         """
        
# rule(s) to download reference data based on the references.tsv file

rule ref_download:
    output:
        ref = "data/ref/{reference_type}/{reference_file}"
    conda:
        "../envs/curl.yaml"
    params:
        file = lambda wildcards: references.loc[wildcards.reference_type].get("file_download"),
    resources:
        parallel_download_connections = 1
    log:
        "logs/ref_download/{reference_type}/{reference_file}.log"
    shell:
        """
        if [[ "{params.file}" =~ \.gz$ ]] && [[ ! "{wildcards.reference_file}" =~ \.gz$ ]]
        then
            ( curl -sS --output {output.ref}.gz {params.file} ) 2> {log}
            ( gzip -d {output.ref}.gz ) 2>> {log}
        else
            ( curl -sS --output {output.ref} {params.file} ) 2> {log}
        fi
        """

ruleorder: create_dict > bwa_index > idx_download > ref_download

rule idx_download:
    output:
        idx = "data/ref/{reference_type}/{reference_file}.{idx_ext}"
    conda:
        "../envs/curl.yaml"
    params:
        index = lambda wildcards: references.loc[wildcards.reference_type].get("index_download")
    resources:
        parallel_download_connections = 1
    log:
        "logs/ref_download/{reference_type}/{reference_file}.{idx_ext}.log"
    shell:
        """
        ( curl -sS --output {output.idx} {params.index} ) 2> {log}
        """


