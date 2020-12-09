import os
import re
import yaml
import argparse
import math
import pandas as pd

workdir: "/fast/users/szyskam_c/scratch/projects/GDCrest/"
bamdir = "/fast/users/szyskam_c/work/NGSData/GDC/orgbam"
sample_sheet = "/fast/users/szyskam_c/snakes/develop/bam2fastq/sheets/GDC_large.txt"

# get the sample_df
sample_df = pd.read_csv(sample_sheet, sep='\t')

sample_list = list(sample_df['sample'])
print(sample_list)


def get_bam(w):
    s = w.sample
    bam = sample_df.set_index('sample').loc[s, 'bam']
    return os.path.join(bamdir, bam)


rule all:
    input: expand("fastq/{sample}_{read}.fastq.gz", read=['R1', 'R2'], sample=sample_list)


# rule bam2fastq:
#     input:
#         bam = get_bam
#     output:
#         fastq1 = "fastq/{sample}_R1.fastq.gz",
#         fastq2 = "fastq/{sample}_R2.fastq.gz"
#     threads:
#         20
#     conda:
#         "align-env.yml"
#     params:
#         deflater = "use_jdk_deflater=true use_jdk_inflater=true"
#     shell:
#         "picard -Xmx60G SamToFastq {params.deflater} "
#         "I={input.bam} F=>(gzip > {output.fastq1}) F2=>(gzip > {output.fastq2}) "
#         "NON_PF=TRUE "
#         "VALIDATION_STRINGENCY=LENIENT "
#         "TMP_DIR=$TMPDIR"


rule bam2fastq:
    input:
        bam = get_bam
    output:
        fastq1 = "fastq/{sample}_R1.fastq.gz",
        fastq2 = "fastq/{sample}_R2.fastq.gz",
        unpaired = "fastq/{sample}_unpaired.fastq.gz"
    threads:
        20
    conda:
        "align-env.yml"
    params:
        deflater = "use_jdk_deflater=true use_jdk_inflater=true"
    shell:
        "samtools sort -n {input.bam} | samtools fastq -0 {output.unpaired} -1 {output.fastq1} -2 {output.fastq2} -@ {threads} -"