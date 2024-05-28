#!/bin/bash

# Job name
#$ -N minimap2 

# Execute the script form the current working directory
#$ -cwd

# Merge the output of the script, and any error messages generated to one file
#$ -j n

# Send the output of the script to a directory called 'UGE-output' in the 
# current working directory (cwd)
# **NOTE: Be sure to create this directory before submitting your job, UGE 
# scheduler will NOT create this direcotry for you**
#$ -o /hpcdata/lvd_qve/Sequencing_Data/QVEU_Seq_0027_Miseq_EVD68_Passage5/

# Tell the job your memory requirements
#$ -l mem_free=10G,h_vmem=12G

# Send mail when the job is submitted, and when the job completes
#$ -m be

# Specify an email address to use
#$ -M adam.catching@nih.gov

#==========================================================================#
# This script is for automating the minimap2 alignment of Nanopore reads
# And for batching the reference sequence by row of a 96-well plate
# Built for the EV-D68 GxE project

module load minimap2/2.10

# Define base directory (now fed in as input arguments)
base_dir=$1

# Define location of reference sequences (now fed in as input arguments)
#fermon_ref="../data/sequences/fermon.fa"
#MO_ref="../data/sequences/MO.fa"
fermon_ref=$2
MO_ref=$3


# Loop through N samples
N=96
# for i in {1..24} --> This is old code for reference
for i in $( eval echo {0..$N} )
do
  
  # This function groups the i-th sample into groups of 8 (for reference sequence)
  alternate=$(( ( (i-1)/8 ) % 2))

  # There is an extra 0 in the name of the first 9 reads that will mess up location names
  if (($i < 10))
    then 
      udp="barcode0"$i
    else 
      udp="barcode"$i
  fi

  # Define the location of the reads for this barcode
  read=$base_dir"barcoding/"$udp"/"$udp".fastq"

  # Define the location of the output files (requires > to be in the string)
  aligned_output=$base_dir"sam/"$udp".sam"

  # Align based on i-th reference sequence
  if (($alternate==0))
      then
        #echo minimap2 -ax map-ont MO.fa $read ">" $aligned_output
        minimap2 -ax map-ont $MO_ref $read > $aligned_output
      else
        #echo minimap2 -ax map-ont fermon.fa $read  ">" $aligned_output
        minimap2 -ax map-ont $fermon_ref $read > $aligned_output
  fi
done

echo All reads have been mapped