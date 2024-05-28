#!/bin/bash

# Job Name
#$ -N Surfseq

#SBATCH --job-name=fastq2gatk
#SBATCH --output=fastq2gatk.out
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64g

#  Specify an email address to use
#$ -M catchingba@nih.gov

#==========================================================================#
# This script is for steps of running Minimap2, convert each .sam file to a 
# Q20 file of aligned basecalls, then 
# EV-D68 genomes. Built for the EV-D68 GxE project

"""
RUN MINIMAP ALONG ALL PASSAGES
"""

# Load modules
module load minimap2

# Change the names of the files relative to where the script is run from 
pass1_dir="/data/home/catchingba/lvd_qve/Sequencing_Data/QVEU_Seq_0046_Minion_BAC_EVD68_ARTIC_passaging_1/no_sample/20230203_2340_MC-113212_FAV75029_16c1fda2/"
pass2_dir="/data/home/catchingba/lvd_qve/Sequencing_Data/QVEU_Seq_0054_Minion_BAC_EVD68_ARTIC_passaging_2/no_sample/20230222_1758_MC-113212_FAV69381_bbde0cc7/"
pass3_dir="/data/home/catchingba/lvd_qve/Sequencing_Data/QVEU_Seq_0059_Minion_BAC_EVD68_ARTIC_passaging_3/no_sample/20230310_2036_MN42546_FAV69394_c5ed8d96/"
pass4_dir="/data/home/catchingba/lvd_qve/Sequencing_Data/QVEU_Seq_0066_Minion_BAC_EVD68_ARTIC_passaging_4/no_sample/20230404_2112_MN42546_FAW92889_67d42911/"
pass5_dir="/data/home/catchingba/lvd_qve/Sequencing_Data/QVEU_Seq_0078_Nanopore_BAC_EVD68_ARTIC_P5/no_sample/20230525_1951_MN42546_FAW93540_f64e53a1/"
pass6_dir="/data/home/catchingba/lvd_qve/Sequencing_Data/QVEU_Seq_0082_Nanopore_BAC_EVD68_ARTIC_P6/no_sample/20230616_1855_MN42546_FAX19209_c0ce6ce4/"
fermon_ref="data/sequences/fermon.fa"
MO_ref="data/sequences/MO.fa"

# Create an array for locations of the passaged sequences
passages=($pass1_dir $pass2_dir $pass3_dir $pass4_dir $pass5_dir $pass6_dir)

# Align several, multiplexed fastq files
echo "Mapping multiplexed reads from all of the reads"
for dir in ${passages[@]}
do
	echo "minimap_96_sample_loop.sh" $dir $fermon_ref $MO_ref
	src/Surfseq/minimap_96_sample_loop.sh $dir $fermon_ref $MO_ref
done

# Done with Minimap, purge
module purge minimap2

"""
RUN SAM2Q20, RESULTING IN BARCODED Q20 FILES IN FOLDER HOLDING ORIGINAL SEQUENCES
"""

module load python
module load numpy

# Define base directory (change as needed)

for dir in ${passages[@]}
do
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
		samfile=$dir"sam/"$udp".sam"
		# Define the location of the output files (requires > to be in the string)
		q20_output=$dir"Q20/"$udp".txt"

		  # Align based on i-th reference sequence
		if (($alternate==0))
			then
				echo "sam2Q20.py" $MO_re f$q20_output
				python src/Surfseq/sam2Q20.py $MO_ref $samfile $q20_output
			else
				echo "sam2Q20.py" $fermon_ref $q20_output
				python src/Surfseq/sam2Q20.py $fermon_ref $samfile $q20_output
			fi
		done

	echo All reads have been mapped

module purge python
module purge numpy

"""
RUN SCULPIN TO RETURN DATAFRAME WITH AMINO ACID SUBSITUTIONS FOR EACH PASSAGE
"""

# Load in modules
module load biopython
module load pandas
module load numpy

python src/Surfseq/sculpin_v0.6.py

# Be kind, rewind
module purge biopython
module purge pandas
module purge numpy
