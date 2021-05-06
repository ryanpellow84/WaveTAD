#!/bin/bash

#### Input Arguements ####
BAM1=$1
BAM2=$2
CHROM_SIZES=$3

#### Create Directory ####
NAM=$(echo $BAM1 | sed 's/\_1.bam//g')

#### Sam Files ####
SAM1=${NAM}_corrected_1.sam
SAM2=${NAM}_corrected_2.sam

#### Bed Files ####
BED1=${NAM}_corrected_1.bed
BED2=${NAM}_corrected_2.bed
BED3=${NAM}_right_contacts.bed
BED4=${NAM}_left_contacts.bed
BED5=${NAM}_right_cov_full.bed
BED6=${NAM}_left_cov_full.bed

#### R Programs ####
R_PROG1=wavelet_splitter.R

#### Coverage Files ####
RIGHT_COV=${NAM}_right_coverage.txt
LEFT_COV=${NAM}_left_coverage.txt

#### Wavelet Splitter ####
samtools view -F 4 $BAM1 | cut -f 1-15 | sed 's/       4       \*      /       4       \*      0       0       /g' | awk '$2 == "0" || $2 == "16"' > $SAM1
cut -f 1-4 $SAM1 > $BED1

samtools view -F 4 $BAM2 | cut -f 1-15 | sed 's/       4       \*      /       4       \*      0       0       /g' | awk '$2 == "0" || $2 == "16"' > $SAM2
cut -f 1-4 $SAM2 > $BED2

rm $SAM2

READ_LEN=$(awk '{print $10}' $SAM1 | head -n 1 | wc -m)

rm $SAM1

R CMD BATCH --no-save --no-restore "--args $BED1 $BED2 $BED3 $BED4 $READ_LEN $CHROM_SIZES" $R_PROG1

rm $BED1
rm $BED2

bedtools genomecov -i $BED3 -g $CHROM_SIZES -d > $BED5
bedtools genomecov -i $BED4 -g $CHROM_SIZES -d > $BED6

rm $BED3
rm $BED4

awk '$3 != "0"' $BED5 > $RIGHT_COV
awk '$3 != "0"' $BED6 > $LEFT_COV

rm $BED5
rm $BED6

echo Calculated Hi-C Coverages




