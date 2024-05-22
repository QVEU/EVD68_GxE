#!/bin/bash

# Job name
#$ -N Sculpin

# Execute the script form the current working directory
#$ -cwd

# Merge the output of the script, and any error messages generated to one file
#$ -j n

# Send the output of the script to a directory called 'UGE-output' in the 
# current working directory (cwd)
# **NOTE: Be sure to create this directory before submitting your job, UGE 
# scheduler will NOT create this direcotry for you**
#$ -o /hpcdata/lvd_qve/Projects/EVD68_adaptation/scripts/sculpin_submit_errors/
#$ -e /hpcdata/lvd_qve/Projects/EVD68_adaptation/scripts/sculpin_submit_errors/
# Tell the job your memory requirements
#$ -l mem_free=12G,h_vmem=12G

# Send mail when the job is submitted, and when the job completes
#$ -m be

# Specify an email address to use
#$ -M adam.catching@naiad.nih.gov

#==========================================================================#
# This script is for automating the generation of a Q20-cutoff of aligned EV-D68 genomes

module load biopython/1.65-foss-2016a-Python-2.7.11
module load pandas/0.18.0-foss-2016a-Python-2.7.11

python scripts/sculpin_v0.4.py 