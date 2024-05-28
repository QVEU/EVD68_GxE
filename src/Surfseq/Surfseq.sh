# Benjamin Adam Catching
# 2023-05-28
# NIH-NIAID-LVD
# Quantitative Virology and Evolution Unit

# Job name
#$ -N Surfseq

# Execute the script form the current working directory
#$ -cwd

# Merge the output of the script, and any error messages generated to one file
#$ -j n

# Send the output of the script to a directory called 'UGE-output' in the 
# current working directory (cwd)
# **NOTE: Be sure to create this directory before submitting your job, UGE 
# scheduler will NOT create this direcotry for you**
#$ -o .
# Tell the job your memory requirements
#$ -l mem_free=10G,h_vmem=12G

# Send mail when the job is submitted, and when the job completes
#$ -m be

# Specify an email address to use
#$ -M adam.catching@naiad.nih.gov

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