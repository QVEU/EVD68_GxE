# Run this for each passage

# Execute from passage folder
mkdir ../sequences/codonvar_per_gene
mkdir ../sequences/dNdS_per_gene

# Define location of reference sequences
fermon_ref="../sequences/fermon.fa"
MO_ref="../sequences/MO.fa"

# Define the location of bed files
fermon_bed="../sequences/fermon_per_gene.bed"
MO_bed="../sequences/MO_per_gene.bed"

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
  read="bam/"$udp".bam"
  # Define codon variation output file
  codonvar_output="codonvar_per_gene/"$udp".csv"
  # Define dnds output file
  dNdS_output="dNdS_per_gene/"$udp".csv"


  # Align based on i-th reference sequence
  if (($alternate==0))
      then
      	echo MO
      	echo Calling $udp codons
        quasitools call codonvar $read $MO_ref 0 $MO_bed -o $codonvar_output
        echo Calling $udp dNdS
        quasitools dnds $codonvar_output $MO_ref 0 > $dNdS_output
      else
      	echo Fermon
      	echo Calling $udp codons
        quasitools call codonvar $read $fermon_ref 0 $fermon_bed -o $codonvar_output
        echo Calling $udp dNdS
        quasitools dnds $codonvar_output $fermon_ref 0 > $dNdS_output
  fi
done

echo All reads have been mapped