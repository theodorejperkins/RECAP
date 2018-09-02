#!/bin/bash
# ===============================================================
# RECAP Wrapper
# RECAP is a wrapper algorithm that resamples ChIP-seq and control
# data to estimate and control for biases built into peak calling
# algorithms. 
# This wrapper script is designed to conveniently peak call 
# ChIP-seq data using MACS and recalibrate peak p-values using RECAP 
# in one convenient package.
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
# Output directory for subsequent MACS and RECAP analyses
OUTPUT_DIR="/global/home/hpc3862/testRECAP"
# Number of remixes for RECAP recalibration (default=1)
BOOTSTRAP=2
# Number of header lines in MACS summary file (default=29 usually)
HEADER=29
# ===============================================================


## PLEASE EDIT AND ADD YOUR DESIRED MACS PARAMETERS IN 2) AND 3)
# ===============================================================
# 1) Re-mix ChIP and control bed files
bash $REMIX_PATH -i $INPUT_DIR -t $CHIP_NAME -c $CONTROL_NAME -o $OUTPUT_DIR -m unequal -b $BOOTSTRAP

# 2) Call original peaks using MACS
# Please specify your own MACS parameters!
# NOTE: p-value threshold must be set to 0.1 for MACS
macs2 callpeak -t $CHIP_NAME -c $CONTROL_NAME --pvalue 0.10 -n ${CHIP_NAME%.*} --outdir "$OUTPUT_DIR/MACS_original"

# 3) Call re-mixed peaks using MACS specifying desired parameters
# Please specify your own MACS parameters!
# NOTE: p-value threshold must be set to 0.1 for MACS
cd "$OUTPUT_DIR/re-mix"
for (( i=1; i<=$BOOTSTRAP; i++ ))
do
	macs2 callpeak -t "${CHIP_NAME%.bed}.bootstrap_$i.bed" -c "${CONTROL_NAME%.bed}.bootstrap_$i.bed" --pvalue 0.10 -n "${CHIP_NAME%.*}.bootstrap_$i" --outdir "$OUTPUT_DIR/MACS_re-mix"
done
# All non-MACS summary files in MACS_re-mix must be deleted if $BOOTSTRAP > 1
cd "$OUTPUT_DIR/MACS_re-mix"
find . -type f ! -name '*_peaks.xls' -delete 

# 4) Recalibrate original peak p-values using RECAP
# NOTE: Check for correct header and p-value column if you obtain any errors here
cd $OUTPUT_DIR
mkdir MACS_RECAP
for (( i=1; i<=$BOOTSTRAP; i++ ))
do
	perl $PERL_PATH --dirOrig "$OUTPUT_DIR/MACS_original" --nameOrig "${CHIP_NAME%.*}_peaks.xls" --dirRemix "$OUTPUT_DIR/MACS_re-mix" --nameRemix "${CHIP_NAME%.*}" --dirOutput "$OUTPUT_DIR/MACS_RECAP" --nameOutput "${CHIP_NAME%.*}.RECAP.bootstrap_${i}_peaks.xls" --bootstrap $i --header $HEADER --pvalCol 7 --delim t --software M
done
# ===============================================================
