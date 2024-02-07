#!/bin/bash
# ===============================================================
# RECAP Wrapper
# RECAP is a wrapper algorithm that resamples ChIP-seq and control
# data to estimate and control for biases built into peak calling
# algorithms. 
# This wrapper script is designed to conveniently peak call 
# ChIP-seq data using diffReps and recalibrate peak p-values using RECAP 
# in one convenient package.
# A detailed explanation of the algorithm can be found under 
# PUBLICATIONS. 
#
# HISTORY:
#   29/08/2017 - v1.0.0 - First Creation
#   14/01/2019 - v1.0.1 - Proper command line parameters
#   15/01/2019 - v1.0.2 - Input validation
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
# Script version number
VERSION="1.0.2"     
# Provide a variable for the location of this and other scripts
SCRIPT_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REMIX_PATH=${SCRIPT_PATH}/RECAP_Re-Mix.sh
PERL_PATH=${SCRIPT_PATH}/RECAP.pl
# Text display commands
bold=$(tput bold)
normal=$(tput sgr0)
# ===============================================================


function mainWrapper() {
####################### Begin Script Here #######################
#################################################################

## PLEASE EDIT AND ADD YOUR DESIRED diffReps PARAMETERS IN 2) AND 3)
# ===============================================================

if [[ ! -d $INPUT_DIR ]]
then
  echo -e "\nERROR: Input directory does not exist"
  exit 1
fi

if [[ ! -d $OUTPUT_DIR ]]
then
  echo -e "\nERROR: Output directory does not exist"
  exit 1
fi

cd $INPUT_DIR

if [[ ! -e $CHIP_NAME  ]]
then
  echo -e "\nERROR: Treatment bed file does not exist"
  exit 1
fi

if [[ ! -e $CONTROL_NAME  ]]
then
  echo -e "\nERROR: Control bed file does not exist."
  echo -e "       Must include the full extension."
  exit 1
fi

if [[ ! $BOOTSTRAP -gt 0 ]]
then
  echo -e "\nERROR: Specify the number of re-mixes"
  echo -e "       Must be a natural number"
  exit 1
fi

if [[ ! $HEADER -ge 0 ]]
then
  echo -e "\nERROR: Specify the number of header lines in output"
  echo -e "       Must be an integer (zero or greater)"
  exit 1
fi



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
  perl $PERL_PATH --dirOrig "$OUTPUT_DIR/diffReps_original" --nameOrig "${CHIP_NAME%.bed}.txt" --dirRemix "$OUTPUT_DIR/diffReps_re-mix" --nameRemix "${CHIP_NAME%.*}" --dirOutput "$OUTPUT_DIR/diffReps_RECAP" --nameOutput "${CHIP_NAME%.bed}.RECAP.bootstrap_${BOOTSTRAP}.txt" --bootstrap $BOOTSTRAP --header $HEADER --pvalCol 13 --delim t --software D


#################################################################
######################## End Script Here ########################
}


#################### Begin Options and Usage ####################
# Print usage
usage() {
  echo -n "
  [Input directory]  [Treatment file] [Control file] 
  [Output directory] [Bootstrap]  [Header]

 ${bold}USAGE:${normal}
  -i, --input 	    Input file directory (absolute path)
  -t, --treatment   Treatment file (full name with extension)
  -c, --control     Control file (full name with extension)
  -o, --output      Output file directory (absolute path)
  -b, --bootstrap   Number of re-mixes
  -e, --header      Header number of peak calling output files

  ${bold}OPTIONS:${normal}
  -h, --help        Display this help and exit
"
}

# ===============================================================
# Iterate over options breaking --foo=bar into --foo bar
unset options
while (($#)); do
  case $1 in
    # If option is of type --foo=bar
    --?*=*) options+=("${1%%=*}" "${1#*=}") ;;
    # add --endopts for --
    --) options+=(--endopts) ;;
    # Otherwise, nothing special
    *) options+=("$1") ;;
  esac
  shift
done
set -- "${options[@]}"
unset options
# ===============================================================


# ===============================================================
# Print help if no arguments or the incorrect number were passed.
[[ $# -eq 0 ]] && set -- "--help"
[[ $# -lt 6 ]] && set -- "--help"

# ===============================================================


# ===============================================================
# Read the options and set stuff
while [[ $1 = -?* ]]; do
  case $1 in
   	-i|--input)     shift; INPUT_DIR=${1} ;;
   	-t|--treatment) shift; CHIP_NAME=${1} ;;
   	-c|--control)   shift; CONTROL_NAME=${1} ;;
   	-o|--output)    shift; OUTPUT_DIR=${1} ;;
	-b|--bootstrap) shift; BOOTSTRAP=${1} ;;
	-e|--header)    shift; HEADER=${1} ;;
	-h|--help)      usage >&2; exit 0 ;;
		*)       echo "ERROR: Bad argument ${1}" ; exit 1 ;;
	esac
	shift
done	

# Store the remaining part as arguments.
args+=("$@")
# ===============================================================


##################### End Options and Usage #####################
# ===============================================================
# Set IFS to preferred implementation
IFS=$'\n\t'
# ===============================================================


# ===============================================================
# Run script
mainWrapper
# ===============================================================
