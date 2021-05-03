#!/bin/bash

#### Input Arguements ####
GZFASTQ1=$1
GZFASTQ2=$2
REF=$3
CHROM_SIZES=$4
OUT_DIR=$5

#### Create Directory ####
NAM=$(echo $GZFASTQ1 | sed 's/\_1.fastq.gz//g')
echo $NAM
BIG_CHR=$(head -n 1 ${CHROM_SIZES} | cut -f 2)
QCFOLDER=$OUT_DIR/qcfolder
mkdir -p QCFOLDER

#### Resolution Name ####
RES_NAM1=1kb
RES_NAM2=5kb
RES_NAM3=10kb
RES_NAM4=25kb
RES_NAM5=50kb

#### FASTQ Files ####
FASTQ1=${NAM}_1.fastq
FASTQ2=${NAM}_2.fastq

#### Bam Files ####
BAM1=${NAM}_1.bam
BAM2=${NAM}_2.bam
BAM3=${NAM}_right_contacts.bam
BAM4=${NAM}_left_contacts.bam

#### Sam Files ####
SAM1=${NAM}_corrected_1.sam
SAM2=${NAM}_corrected_2.sam
SAM3=${NAM}_right_contacts_no_header.sam
SAM4=${NAM}_left_contacts_no_header.sam
SAM5=${NAM}_right_contacts.sam
SAM6=${NAM}_left_contacts.sam

#### Bed Files ####
BED1=${NAM}_corrected_1.bed
BED2=${NAM}_corrected_2.bed
BED3=${NAM}_right_contacts.bed
BED4=${NAM}_left_contacts.bed
BED5=${NAM}_right_cov_full.bed
BED6=${NAM}_left_cov_full.bed

#### R Programs ####
R_PROG1=wavelet_splitter.R

#### Cool Files ####
COOL1=${NAM}_1kb.cool
COOL2=${NAM}_5kb.cool
COOL3=${NAM}_10kb.cool
COOL4=${NAM}_25kb.cool
COOL5=${NAM}_50kb.cool
CORRECT1=${NAM}_1kb_corrected.cool
CORRECT2=${NAM}_5kb_corrected.cool
CORRECT3=${NAM}_10kb_corrected.cool
CORRECT4=${NAM}_25kb_corrected.cool
CORRECT5=${NAM}_50kb_corrected.cool
NORM1=${NAM}_1kb_norm_corrected.cool
NORM2=${NAM}_5kb_norm_corrected.cool
NORM3=${NAM}_10kb_norm_corrected.cool
NORM4=${NAM}_25kb_norm_corrected.cool
NORM5=${NAM}_50kb_norm_corrected.cool

#### Contact Files ####
CONTACTS1=${OUT_DIR}/${NAM}_1kb_contacts.txt
CONTACTS2=${OUT_DIR}/${NAM}_5kb_contacts.txt
CONTACTS3=${OUT_DIR}/${NAM}_10kb_contacts.txt
CONTACTS4=${OUT_DIR}/${NAM}_25kb_contacts.txt
CONTACTS5=${OUT_DIR}/${NAM}_50kb_contacts.txt

#### Miscellaneous Files ####
HEADER=${OUT_DIR}/${NAM}_header.txt
RIGHT_COV=${NAM}_right_coverage.txt
LEFT_COV=${NAM}_left_coverage.txt
MAD_VAL1=${NAM}_1kb_mad_threshold_value.txt
MAD_VAL2=${NAM}_5kb_mad_threshold_value.txt
MAD_VAL3=${NAM}_10kb_mad_threshold_value.txt
MAD_VAL4=${NAM}_25kb_mad_threshold_value.txt
MAD_VAL5=${NAM}_50kb_mad_threshold_value.txt
CHR_ONLY=${NAM}_chromosomes_only.txt
awk '{print $1}' $CHROM_SIZES > $CHR_ONLY
readarray CHR_ARR < $CHR_ONLY

#### Bin Sizes ####
BIN_SIZE1=1000
BIN_SIZE2=5000
BIN_SIZE3=10000
BIN_SIZE4=25000
BIN_SIZE5=50000

#### Build Matrix Logs ####
BUILD_LOG1=${NAM}_1kb_build_mat_qc.log
BUILD_LOG2=${NAM}_5kb_build_mat_qc.log
BUILD_LOG3=${NAM}_10kb_build_mat_qc.log
BUILD_LOG4=${NAM}_25kb_build_mat_qc.log
BUILD_LOG5=${NAM}_50kb_build_mat_qc.log

#### Hi-C Explorer Loops ####
LOOP_COOL1=${COOL1}
LOOP_COOL2=${COOL2}
LOOP_COOL3=${COOL3}
LOOP_COOL4=${COOL4}
LOOP_COOL5=${COOL5}
LOOPS1=${OUT_DIR}/${NAM}_1kb_loops_hicexplorer.bedgraph
LOOPS2=${OUT_DIR}/${NAM}_5kb_loops_hicexplorer.bedgraph
LOOPS3=${OUT_DIR}/${NAM}_10kb_loops_hicexplorer.bedgraph
LOOPS4=${OUT_DIR}/${NAM}_25kb_loops_hicexplorer.bedgraph
LOOPS5=${OUT_DIR}/${NAM}_50kb_loops_hicexplorer.bedgraph

#### TopDom Domains ####
TEMP_MAT1=${OUT_DIR}/${NAM}_1kb_temp_hic_matrix.bed
TEMP_MAT2=${OUT_DIR}/${NAM}_5kb_temp_hic_matrix.bed
TEMP_MAT3=${OUT_DIR}/${NAM}_10kb_temp_hic_matrix.bed
TEMP_MAT4=${OUT_DIR}/${NAM}_25kb_temp_hic_matrix.bed
TEMP_MAT5=${OUT_DIR}/${NAM}_50kb_temp_hic_matrix.bed
TOPDOM1=${OUT_DIR}/${NAM}_1kb_tad_domains_topdom.bed
TOPDOM2=${OUT_DIR}/${NAM}_5kb_tad_domains_topdom.bed
TOPDOM3=${OUT_DIR}/${NAM}_10kb_tad_domains_topdom.bed
TOPDOM4=${OUT_DIR}/${NAM}_25kb_tad_domains_topdom.bed
TOPDOM5=${OUT_DIR}/${NAM}_50kb_tad_domains_topdom.bed

#### WUBBA Domains ####
WUBBA=${OUT_DIR}/${NAM}_TAD_domains_wubba_topdom_collab_results.bed

#### BWA ####
####gunzip $GZFASTQ1
####gunzip $GZFASTQ2

####echo $REF
####echo $FASTQ1
####echo $BAM1

####bwa mem -t 56 -E 50 -L 0 $REF $FASTQ1 | samtools view --threads 56 -bS - -o $BAM1
####bwa mem -t 56 -E 50 -L 0 $REF $FASTQ2 | samtools view --threads 56 -bS - -o $BAM2

####samtools view $BAM1 > $OUT_DIR/sam_test.sam

#### Wavelet Splitter ####
####samtools view -F 4 $BAM1 | cut -f 1-15 | sed 's/       4       \*      /       4       \*      0       0       /g' | awk '$2 == "0" || $2 == "16"' > $SAM1
####cut -f 1-4 $SAM1 > $BED1

####samtools view -F 4 $BAM2 | cut -f 1-15 | sed 's/       4       \*      /       4       \*      0       0       /g' | awk '$2 == "0" || $2 == "16"' > $SAM2
####cut -f 1-4 $SAM2 > $BED2

####READ_LEN=$(awk '{print $10}' $SAM1 | head -n 1 | wc -m)

####R CMD BATCH --no-save --no-restore "--args $BED1 $BED2 $BED3 $BED4 $READ_LEN $CHROM_SIZES" $R_PROG1

####bedtools genomecov -i $BED3 -g $CHROM_SIZES -d > $BED5
####bedtools genomecov -i $BED4 -g $CHROM_SIZES -d > $BED6

####awk '$3 != "0"' $BED5 > $RIGHT_COV
####awk '$3 != "0"' $BED6 > $LEFT_COV

#### Build Matrix ####
####if (($BIG_CHR < 35000000)); then
####	hicBuildMatrix -s $BAM1 $BAM2 -bs $BIN_SIZE1 -o $COOL1 --skipDuplicationCheck --QCfolder $QCFOLDER --threads 8
####fi
####hicBuildMatrix -s $BAM1 $BAM2 -bs $BIN_SIZE2 -o $COOL2 --skipDuplicationCheck --QCfolder $QCFOLDER --threads 8
####hicBuildMatrix -s $BAM1 $BAM2 -bs $BIN_SIZE3 -o $COOL3 --skipDuplicationCheck --QCfolder $QCFOLDER --threads 8
####hicBuildMatrix -s $BAM1 $BAM2 -bs $BIN_SIZE4 -o $COOL4 --skipDuplicationCheck --QCfolder $QCFOLDER --threads 8
hicBuildMatrix -s $BAM1 $BAM2 -bs $BIN_SIZE5 -o $COOL5 --skipDuplicationCheck --QCfolder $QCFOLDER --threads 8 -r 3R:25000000-31000000
	
#### Correct Matrix ####
####readarray CHR_ARR < $CHR_ONLY
####if (($BIG_CHR < 35000000)); then
####	R CMD BATCH --no-save --no-restore "--args $BUILD_LOG1 $MAD_VAL1" $R_PROG
####	VAL1=$(head -n 1 $MAD_VAL1)
####	hicCorrectMatrix correct --matrix $OUT_DIR/$COOL1 --filterThreshold $VAL1 5 --chromosomes ${CHR_ARR[@]} -o $CORRECT1
####fi
####R CMD BATCH --no-save --no-restore "--args $BUILD_LOG2 $MAD_VAL2" $R_PROG
####R CMD BATCH --no-save --no-restore "--args $BUILD_LOG3 $MAD_VAL3" $R_PROG
####R CMD BATCH --no-save --no-restore "--args $BUILD_LOG4 $MAD_VAL4" $R_PROG
####R CMD BATCH --no-save --no-restore "--args $BUILD_LOG5 $MAD_VAL5" $R_PROG
####VAL2=$(head -n 1 $MAD_VAL2)
####VAL3=$(head -n 1 $MAD_VAL3)
####VAL4=$(head -n 1 $MAD_VAL4)
####VAL5=$(head -n 1 $MAD_VAL5)
####hicCorrectMatrix correct --matrix $OUT_DIR/$COOL2 --filterThreshold $VAL1 5 --chromosomes ${CHR_ARR[@]} -o $CORRECT2
####hicCorrectMatrix correct --matrix $OUT_DIR/$COOL3 --filterThreshold $VAL1 5 --chromosomes ${CHR_ARR[@]} -o $CORRECT3
####hicCorrectMatrix correct --matrix $OUT_DIR/$COOL4 --filterThreshold $VAL1 5 --chromosomes ${CHR_ARR[@]} -o $CORRECT4
####hicCorrectMatrix correct --matrix $OUT_DIR/$COOL5 --filterThreshold $VAL1 5 --chromosomes ${CHR_ARR[@]} -o $CORRECT5

#### Normalize Matrix ####
####if (($BIG_CHR < 35000000)); then
####	hicNormalize -m $CORRECT1 --normalize norm_range -o $NORM1
####fi
####hicNormalize -m $CORRECT2 --normalize norm_range -o $NORM2
####hicNormalize -m $CORRECT3 --normalize norm_range -o $NORM3
####hicNormalize -m $CORRECT4 --normalize norm_range -o $NORM4
####hicNormalize -m $CORRECT5 --normalize norm_range -o $NORM5

#### Detect Loops ####
####if (($BIG_CHR < 35000000)); then
####hicDetectLoops -m $LOOP_COOL1 -o $LOOPS1 --chromosomes ${CHR_ARR[@]} -p 1 -pw 2 -w 5 -pp 0.1 -pit 10 -oet 1.5 --maxLoopDistance 5000000
####fi
####hicDetectLoops -m $LOOP_COOL2 -o $LOOPS2 --chromosomes ${CHR_ARR[@]} -p 1 -pw 2 -w 5 -pp 0.1 -pit 10 -oet 1.5 --maxLoopDistance 5000000
####hicDetectLoops -m $LOOP_COOL3 -o $LOOPS3 --chromosomes ${CHR_ARR[@]} -p 1 -pw 2 -w 5 -pp 0.1 -pit 10 -oet 1.5 --maxLoopDistance 5000000
####hicDetectLoops -m $LOOP_COOL4 -o $LOOPS4 --chromosomes ${CHR_ARR[@]} -p 1 -pw 2 -w 5 -pp 0.1 -pit 10 -oet 1.5 --maxLoopDistance 5000000
####hicDetectLoops -m $LOOP_COOL5 -o $LOOPS5 --chromosomes ${CHR_ARR[@]} -p 1 -pw 2 -w 5 -pp 0.1 -pit 10 -oet 1.5 --maxLoopDistance 5000000

#### TopDom ####
####if (($BIG_CHR < 35000000)); then
####	R CMD BATCH --no-save --no-restore "--args $CONTACTS1 $TEMP_MAT1 $TOPDOM1 $BIN_SIZE1 $CHROM_SIZES" $R_PROG
####fi
####R CMD BATCH --no-save --no-restore "--args $CONTACTS2 $TEMP_MAT2 $TOPDOM2 $BIN_SIZE2 $CHROM_SIZES" $R_PROG
####R CMD BATCH --no-save --no-restore "--args $CONTACTS3 $TEMP_MAT3 $TOPDOM3 $BIN_SIZE3 $CHROM_SIZES" $R_PROG
####R CMD BATCH --no-save --no-restore "--args $CONTACTS4 $TEMP_MAT4 $TOPDOM4 $BIN_SIZE4 $CHROM_SIZES" $R_PROG
####R CMD BATCH --no-save --no-restore "--args $CONTACTS5 $TEMP_MAT5 $TOPDOM5 $BIN_SIZE5 $CHROM_SIZES" $R_PROG
	
#### WUBBA-TopDom Collab ####

####JOB_ID_VEC=()
####PASS_FILE_VEC=()
####WUBBA_RESULT_VEC=()
####for i in ${CHR_ARR[@]}
####do
####	WUBBA_RESULT=${NAM}_TAD_domains_wubba_topdom_results_${i}.bed
####	WUBBA_RESULT_VEC+=$WUBBA_RESULT
####	if (($BIG_CHR < 35000000)); then
####		GEN_SIZE=small
####		awk -v myvar="${i}" '$1==myvar' $LEFT_COV > $OUT_DIR/${NAM}_${i}_left_coverage.txt
####		awk -v myvar="${i}" '$1==myvar' $RIGHT_COV > $OUT_DIR/${NAM}_${i}_right_coverage.txt
####		R CMD BATCH --no-save --no-restore "--args $OUT_DIR/${NAM}_${i}_left_coverage.txt $OUT_DIR/${NAM}_${i}_right_coverage.txt $WUBBA_RESULT $TOPDOM1 $TOPDOM2 $TOPDOM3 $TOPDOM4 $LOOPS1 $LOOPS2 $LOOPS3 $LOOPS4 ${i} $GEN_SIZE" $R_PROG ${NAM}_wubba_topdom_collab_${i}.Rout

####	else
####		GEN_SIZE=big
####		R CMD BATCH --no-save --no-restore "--args $OUT_DIR/${NAM}_${i}_left_coverage.txt $OUT_DIR/${NAM}_${i}_right_coverage.txt $WUBBA_RESULT $TOPDOM2 $TOPDOM3 $TOPDOM4 $TOPDOM5 $LOOPS2 $LOOPS3 $LOOPS4 $LOOPS5 ${i} $GEN_SIZE" $R_PROG ${NAM}_wubba_topdom_collab_${i}.Rout
####	fi
####done
	
#### Merge WaveTAD Results ####
####R CMD BATCH --no-save --no-restore "--args $BEST_FILES $BEST" $R_PROG
	





