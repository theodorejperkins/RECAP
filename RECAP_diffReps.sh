#!/bin/bash
# ===============================================================
# RECAP Wrapper
# RECAP is a wrapper algorithm that resamples ChIP-seq and control
# data to estimate and control for biases built into peak calling
# algorithms. 
# This wrapper script is designed to conveniently peak call 
# ChIP-seq data using diffReps and recalibrate peak p-values using 
# RECAP  in one convenient package.
# A detailed explanation of the algorithm can be found under 
# PUBLICATIONS. 
#
# HISTORY:
#   29/08/2017 - v1.0.0 - First Creation
#
# CREDITS:
# RECAP was developed by Justin G. Chitpin, Aseel Awdeh, 
# and Theodore J. Perkins.
# Development of RECAP was carried out at the Ottawa Hospital
# Research Institute in the Perkins Lab.
#
# PUBLICATIONS:
# If you use RECAP, please cite the following paper:
# <INSERT PUBLICATION HERE>
#
# QUESTIONS:
# Please contact tperkins@ohri.ca
# ===============================================================


# ===============================================================
# Provide a variable for the location of this and other scripts
SCRIPT_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REMIX_PATH=$(find ~/ -type f -name "RECAP_Re-Mix.sh" | head -n 1)
PERL_PATH=$(find ~/ -type f -name "RECAP.pl" | head -n 1)
# ===============================================================


## PLEASE FILL THE FOLLOWING PARAMETERS
# ===============================================================
# ChIP/Control directory
INPUT_DIR="/global/home/hpc3862/testRECAP"
# ChIP name
CHIP_NAME="a.bed"
# Control name
CONTROL_NAME="b.bed"
# Output directory for subsequent diffReps and RECAP analyses
OUTPUT_DIR="/global/home/hpc3862/testRECAP"
# Number of remixes for RECAP recalibration (default=1)
BOOTSTRAP=1
# Number of header lines in diffReps summary file (default=33)
HEADER=33
# ===============================================================


## PLEASE EDIT AND ADD YOUR DESIRED diffReps PARAMETERS IN 2) AND 3)
# ===============================================================
# 1) Re-mix ChIP and control bed files
bash $REMIX_PATH -i $INPUT_DIR -t $CHIP_NAME -c $CONTROL_NAME -o $OUTPUT_DIR -m unequal -b $BOOTSTRAP

# 2) Call original peaks using diffReps
# Please specify your own diffReps parameters!
# NOTE: p-value threshold must be set to 1.0 for diffReps
cd $OUTPUT_DIR
mkdir -p diffReps_original
diffReps.pl --treatment "$INPUT_DIR/$CHIP_NAME" --control "$INPUT_DIR/$CONTROL_NAME" -meth gt -gname hg19 --pval 1 --mode n --report "$OUTPUT_DIR/diffReps_original/${CHIP_NAME%.bed}.txt" --nohs --noanno --nsd 20

# 3) Call re-mixed peaks using diffReps specifying desired parameters
# Please specify your own diffReps parameters!
# NOTE: p-value threshold must be set to 1.0 for diffReps
cd $OUTPUT_DIR
mkdir -p diffReps_re-mix
for (( i=1; i<=$BOOTSTRAP; i++ ))
do
	diffReps.pl --treatment "$OUTPUT_DIR/re-mix/${CHIP_NAME%.bed}.bootstrap_$i.bed" --control "$OUTPUT_DIR/re-mix/${CONTROL_NAME%.bed}.bootstrap_$i.bed" -meth gt -gname hg19 --pval 1 --mode n --report "$OUTPUT_DIR/diffReps_re-mix/${CHIP_NAME%.bed}.bootstrap_$i.txt" --nohs --noanno --nsd 20
done

# 4) Recalibrate original peak p-values using RECAP
# NOTE: Check for correct header and p-value column if you obtain any errors here
cd $OUTPUT_DIR
mkdir diffReps_RECAP
for (( i=1; i<=$BOOTSTRAP; i++ ))
do
	perl $PERL_PATH --dirOrig "$OUTPUT_DIR/diffReps_original" --nameOrig "${CHIP_NAME%.bed}.txt" --dirRemix "$OUTPUT_DIR/diffReps_re-mix" --nameRemix "${CHIP_NAME%.bed}.bootstrap_$i.txt" --dirOutput "$OUTPUT_DIR/diffReps_RECAP" --nameOutput "${CHIP_NAME%.bed}.RECAP.bootstrap_$i.txt" --bootstrap $i --header $HEADER --pvalCol 13 --delim t --software D
done
# ===============================================================
