#!/bin/bash

#$ -V
#$ -j y
#$ -N FCAHsub
#$ -o logs
#$ -r yes
#$ -cwd
#$ -S /bin/bash
#$ -P control

export LOGDIR=/fast/users/${USER}/scratch/lox/Bam2Fastq/${JOB_ID}
export TMPDIR=/fast/users/${USER}/scratch/tmp
# export WRKDIR=$HOME/work/projects/whWES

mkdir -p $LOGDIR

# somehow my environments are not set
# have to set it explicitly
# conda activate somVar-EB
conda activate WES-env
# outputs every output to the terminal
set -x

# !!! leading white space is important
DRMAA=" -pe smp {threads}  -l h_rt=04:00:00 -l h_vmem=3.5g"
DRMAA="$DRMAA -V -o $LOGDIR/ -j yes"
snakemake --unlock
snakemake --use-conda --rerun-incomplete --drmaa "$DRMAA" -j 2000 -p -r -k
# -k ..keep going if job fails
# -p ..print out shell commands
