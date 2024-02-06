#!/bin/bash
# ===============================================================
# RECAP Wrapper
# RECAP is a wrapper algorithm that resamples ChIP-seq and control
# data to estimate and control for biases built into peak calling
# algorithms. 
# This wrapper script is designed to conveniently peak call 
# ChIP-seq data using SICER and recalibrate peak p-values using RECAP 
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
REMIX_PATH=$(find ~/ -type f -name "RECAP_Re-Mix.sh" | head -n 1)
PERL_PATH=$(find ~/ -type f -name "RECAP.pl" | head -n 1)
SICER_PATH=$(find ~/ -type f -name "SICER.sh" | head -n 1)
# Text display commands
bold=$(tput bold)
normal=$(tput sgr0)
# ===============================================================


function mainWrapper() {
####################### Begin Script Here #######################
#################################################################

## PLEASE EDIT AND ADD YOUR DESIRED SICER PARAMETERS IN 2) AND 3)
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

# 2) Call original peaks using SICER
# Please specify your own SICER parameters!
# NOTE: p-value threshold must be set to 1.0 for SICER
cd $OUTPUT_DIR
mkdir -p SICER_original
cd $INPUT_DIR
bash $SICER_PATH  $INPUT_DIR $CHIP_NAME $CONTROL_NAME "$OUTPUT_DIR/SICER_original" hg38 1 $WINDOW 100 1 $GAP 1

# 3) Call re-mixed peaks using SICER specifying desired parameters
# Please specify your own SICER parameters!
# NOTE: p-value threshold must be set to 1.0 for SICER
cd $OUTPUT_DIR
mkdir -p SICER_re-mix
cd SICER_re-mix
for (( i=1; i<=$BOOTSTRAP; i++ ))
do
  bash $SICER_PATH  "$OUTPUT_DIR/re-mix" "${CHIP_NAME%.bed}.bootstrap_$i.bed" "${CONTROL_NAME%.bed}.bootstrap_$i.bed" "$OUTPUT_DIR/SICER_re-mix" hg38 1 $WINDOW 100 1 $GAP 1
done

# All non-SICER summary files in SICER_re-mix must be deleted if $BOOTSTRAP > 1
if [ -d "$OUTPUT_DIR/SICER_re-mix" ]
then 
  cd "$OUTPUT_DIR/SICER_re-mix"
  find . -type f ! -name '*-islands-summary' -delete
else
  echo "Output directory doesn't exist!"
  exit
fi

# 4) Recalibrate original peak p-values using RECAP
# NOTE: Check for correct header and p-value column if you obtain any errors here
cd $OUTPUT_DIR
mkdir -p SICER_RECAP

  perl $PERL_PATH --dirOrig "$OUTPUT_DIR/SICER_original" --nameOrig "${CHIP_NAME%.bed}-W${WINDOW}-G${GAP}-islands-summary" --dirRemix "$OUTPUT_DIR/SICER_re-mix" --nameRemix "${CHIP_NAME%.bed}" --dirOutput "$OUTPUT_DIR/SICER_RECAP" --nameOutput "${CHIP_NAME%.bed}.RECAP.bootstrap_${BOOTSTRAP}-W${WINDOW}-G${GAP}-islands-summary" --bootstrap $BOOTSTRAP --header $HEADER --pvalCol 6 --delim t --software O

#################################################################
######################## End Script Here ########################
}


#################### Begin Options and Usage ####################
# Print usage
usage() {
  echo -n "
  [Input directory]  [Treatment file] [Control file] 
  [Output directory] [Window] [Gap] [Bootstrap]  [Header]

 ${bold}USAGE:${normal}
  -i, --input 	    Input file directory (absolute path)
  -t, --treatment   Treatment file (full name with extension)
  -c, --control     Control file (full name with extension)
  -o, --output      Output file directory (absolute path)
  -w, --window      SICER Window size
  -g, --gap         SICER Gap size
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
[[ $# -lt 8 ]] && set -- "--help"

# ===============================================================

# ===============================================================
# Read the options and set stuff
while [[ $1 = -?* ]]; do
  case $1 in
   	-i|--input)     shift; INPUT_DIR=${1} ;;
   	-t|--treatment) shift; CHIP_NAME=${1} ;;
   	-c|--control)   shift; CONTROL_NAME=${1} ;;
   	-o|--output)    shift; OUTPUT_DIR=${1} ;;
	-w|--window)    shift; WINDOW=${1} ;;
	-g|--gap)       shift; GAP=${1} ;;
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
