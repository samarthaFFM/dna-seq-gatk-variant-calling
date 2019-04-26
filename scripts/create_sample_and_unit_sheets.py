import enasearch
'''
The above import requires enasearch to be installed, e.g. using bioconda:
conda create -n enasearch enasearch
conda activate enasearch
'''
import csv

accession = "SRP114315"

wanted_fields = [
        "sample_alias",
        "run_accession",
        "instrument_platform",
        "library_strategy",
        "fastq_ftp",
        "fastq_md5",
        ]

output_fields = [
        "sample",
        "unit",
        "platform",
        "library_strategy",
        "fq1",
        "fq2",
        "fastq_ftp_1",
        "fastq_ftp_2",
        "fastq_md5_1",
        "fastq_md5_2"
        ]
        
samples_file = '../samples.tsv'
units_file = '../units.tsv'

results_dict = csv.DictReader(
                enasearch.retrieve_run_report(
                    accession = accession,
                    fields = ",".join( wanted_fields )
                    ).split('\n'),
                delimiter='\t')

with open(samples_file, 'w') as samples, open(units_file, 'w') as units:
        sample_writer = csv.writer(samples, delimiter = "\t")
        units_writer = csv.writer(units, delimiter = "\t")
        sample_writer.writerow(["sample", "sample_ID", "library_strategy" ])
        units_writer.writerow( output_fields )
        for row in results_dict:
            if row['library_strategy'] == 'WXS':
                sample_ID = row['sample_alias']
                sample_writer.writerow( [ row['sample_alias'], sample_ID, row['library_strategy'] ] )
                (fastq_ftp_1, fastq_ftp_2) = row['fastq_ftp'].split(';')
                (fastq_md5_1, fastq_md5_2) = row['fastq_ftp'].split(';')
                units_writer.writerow(
                    [   row['sample_alias'],
                        row['run_accession'],
                        row['instrument_platform'],
                        row['library_strategy'],
                        'data/raw/' + row['sample_alias'] + '.1.fq.gz',
                        'data/raw/' + row['sample_alias'] + '.2.fq.gz',
                        fastq_ftp_1,
                        fastq_ftp_2,
                        fastq_md5_1,
                        fastq_md5_2
                    ] )

