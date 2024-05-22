#!/bin/bash

# This script performs the Kallisto pipeline on the input fastq files

# Job Name
#$ -N kallisto_index

# Execute the script from the Current Working Directory
#$ -cwd

# Merge the output of the script, and any error messages generated to one file
#$ -j y

# Send the output of the script to a directory called 'UGE-output' in the current working directory (cwd)
# **NOTE: Be sure to create this directory before submitting your job, UGE scheduler will NOT create this**
# **directory for you**
#$ -o /hpcdata/lvd_qve/Sequencing_Data/QVEU_Seq_0109_NextSeq_BAC_Fermon_RD_A549/231208_VH01023_36_AAC2LGHM5
#$ -e /hpcdata/lvd_qve/Sequencing_Data/QVEU_Seq_0109_NextSeq_BAC_Fermon_RD_A549/231208_VH01023_36_AAC2LGHM5

# Tell the job your memory requirements
#$ -l h_vmem=240G

# Send mail when the job is submitted, and when the job completes
#$ -m be

#  Specify an email address to use
#$ -M catchingba@nih.gov

# Download reference genome from https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna
# Choose file Homo_sapiens.GRCh38.cdna.all.fa.gz
module load kallisto
kallisto index Homo_sapiens.GRCh38.cdna.all.fa.gz -i human.idx