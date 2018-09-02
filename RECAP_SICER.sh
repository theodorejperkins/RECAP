#!/bin/bash
# ===============================================================
# RECAP Wrapper
# RECAP is a wrapper algorithm that resamples ChIP-seq and control
# data to estimate and control for biases built into peak calling
# algorithms. 
# This wrapper script is designed to conveniently peak call 
# ChIP-seq data using SICER and recalibrate peak p-values using 
# RECAP in one convenient package.
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
SICER_PATH=$(find ~/ -type f -name "SICER.sh" | head -n 1)
# ===============================================================


## PLEASE FILL THE FOLLOWING PARAMETERS
# ===============================================================
# ChIP/Control directory
INPUT_DIR="/global/home/hpc3862/testRECAP"
# ChIP name
CHIP_NAME="a.bed"
# Control name
CONTROL_NAME="b.bed"
# Output directory for subsequent SICER and RECAP analyses
OUTPUT_DIR="/global/home/hpc3862/testRECAP"
# Number of remixes for RECAP recalibration (default=1)
BOOTSTRAP=2
# Number of header lines in SICER summary file (default=0)
HEADER=0
# ===============================================================


## PLEASE EDIT AND ADD YOUR DESIRED SICER PARAMETERS IN 2) AND 3)
# ===============================================================
# 1) Re-mix ChIP and control bed files
bash $REMIX_PATH -i $INPUT_DIR -t $CHIP_NAME -c $CONTROL_NAME -o $OUTPUT_DIR -m unequal -b $BOOTSTRAP

# 2) Call original peaks using SICER
# Please specify your own SICER parameters!
# NOTE: p-value threshold must be set to 1.0 for SICER
bash $SICER_PATH  $INPUT_DIR $CHIP_NAME $CONTROL_NAME "$OUTPUT_DIR/SICER_original" hg38 1 100 100 1 100 1

# 3) Call re-mixed peaks using SICER specifying desired parameters
# Please specify your own SICER parameters!
# NOTE: p-value threshold must be set to 1.0 for SICER
cd "$OUTPUT_DIR/re-mix"
for (( i=1; i<=$BOOTSTRAP; i++ ))
do
	bash $SICER_PATH  "$OUTPUT_DIR/SICER_re-mix" "${CHIP_NAME%.bed}.bootstrap_$i.bed" "${CONTROL_NAME%.bed}.bootstrap_$i.bed" "$OUTPUT_DIR/SICER_original" hg38 1 100 100 1 100 1
done

# All non-SICER summary files in SICER_re-mix must be deleted if $BOOTSTRAP > 1
cd "$OUTPUT_DIR/SICER_re-mix"
find . -type f ! -name '*-islands-summary' -delete 

# 4) Recalibrate original peak p-values using RECAP
# NOTE: Check for correct header and p-value column if you obtain any errors here
cd $OUTPUT_DIR
mkdir SICER_RECAP
for (( i=1; i<=$BOOTSTRAP; i++ ))
do
	perl $PERL_PATH --dirOrig "$OUTPUT_DIR/SICER_original" --nameOrig "${CHIP_NAME}*" --dirRemix "$OUTPUT_DIR/SICER_re-mix" --nameRemix $CHIP_NAME --dirOutput "$OUTPUT_DIR/SICER_RECAP" --nameOutput "${CHIP_NAME}.RECAP.bootstrap_${i}" --bootstrap $i --header $HEADER --pvalCol 6 --delim t --software O
done
# ===============================================================
