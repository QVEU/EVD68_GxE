# Benjamin Adam Catching
# 2024-05-23
# NIH-NIAID-LVD
# Quantitative Virology and Evolution Unit

#!/bin/bash
#SBATCH --job-name=bam_align
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
#SBATCH -e map.out
module load samtools

# Make the directory of the bam files
# Define base directory

mkdir bam

# Define location of reference sequences
fermon_ref="../sequences/fermon.fa"
MO_ref="../sequences/MO.fa"
samtools faidx $fermon_ref
samtools faidx $MO_ref

# Loop through N samples
N=96
# for i in {1..24} --> This is old code for reference
for i in $( eval echo {1..$N} )
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
  samfile="sam/"$udp".sam"

  # Define the location of the output files (requires > to be in the string)
  bamfile="bam/"$udp".bam"
  bamalign="bam/"$udp".bai"

  # Align based on i-th reference sequence
  samtools view -bS $samfile | samtools sort -o $bamfile
  samtools index $bamfile $bamalign
done

echo All reads have been mapped