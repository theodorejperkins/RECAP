#!/bin/bash
# ===============================================================
# RECAP Re-Mix
# RECAP is a wrapper algorithm that resamples ChIP-seq and control
# data to estimate and control for biases built into peak calling
# algorithms. 
# The purpose of this script is to randomly re-mix ChIP-seq and 
# control prior to peak calling analysis. A detailed explanation
# of the algorithm can be found under PUBLICATIONS. 
#
# HISTORY:
#   03/11/2017 - v1.0.0 - First Creation
#	29/01/2018 - v1.0.1 - Minor commenting update 
#
# CREDITS:
# RECAP was developed by Justin G. Chitpin and Theodore J. Perkins.
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
VERSION="1.0.1"     
# Provide a variable with the location of this script.
SCRIPT_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Text display commands
bold=$(tput bold)
normal=$(tput sgr0)
# ===============================================================


# ===============================================================
# Set Flags
# Flags which can be overridden by user input.
# Default values are below
# ===============================================================


function mainRemix() {
####################### Begin Script Here #######################
#################################################################

start_time=`date +%s`

echo "##################################################"
echo "######         RECAP RE-MIX  v$VERSION         ######"
echo "##################################################"
echo ""
echo "Input directory:    $INPUT_DIR"
echo "Treatment library:  $TREATMENT_NAME"
echo "Control library:    $CONTROL_NAME"
echo "Output directory:   $OUTPUT_DIR"
echo "Re-mix method:      $METHOD_NAME"
echo "Number of re-mixes: $BOOTSTRAP"

if [[ ! -d $INPUT_DIR ]]
then
	echo -e "\nERROR: Input directory does not exist"
	exit 1
fi

cd $INPUT_DIR

if [[ ! -e $TREATMENT_NAME || ! $TREATMENT_NAME == *.bed ]]
then
	echo -e "\nERROR: Treatment bed file does not exist"
	exit 1
elif [[ ! -e $CONTROL_NAME ||  ! $CONTROL_NAME == *.bed ]]
then
	echo -e "\nERROR: Control bed file does not exist"
	exit 1
elif [[ ! -d $OUTPUT_DIR ]]
then
	echo -e "\nERROR: Output directory does not exist"
	exit 1
elif [[ ! $BOOTSTRAP -gt 0 ]]
then
	echo -e "\nERROR: Specify the number of re-mixes"
	echo -e "       Must be a natural number"
	exit 1
fi

# Base names of treatment and control libraries
TREATMENT_NAME_BASE="${TREATMENT_NAME%.*}"
CONTROL_NAME_BASE="${CONTROL_NAME%.*}"
	
for (( BOOTSTRAP_COUNT=1; BOOTSTRAP_COUNT<=$BOOTSTRAP; BOOTSTRAP_COUNT++ ))
do
	cd $INPUT_DIR
	echo ""
	echo "${bold}Starting Re-Mixing Procedure #$BOOTSTRAP_COUNT ${normal}"
	echo "Creating directory for re-mix files:"
	mkdir -p $OUTPUT_DIR/re-mix
	echo "Check!"
	
	echo "Concatenating treatment and control libraries"
	cat $TREATMENT_NAME $CONTROL_NAME > $OUTPUT_DIR/re-mix/$TREATMENT_NAME_BASE.bootstrap_$BOOTSTRAP_COUNT.tmp
	echo "Check!"
	
	# If method is unequal
	if [ $METHOD = 0 ]
	then
		echo "Unequal mixing method selected"
		FIRST_LINES=$( wc -l < $TREATMENT_NAME )
		LAST_LINES=$( wc -l < $CONTROL_NAME )
		cd $OUTPUT_DIR/re-mix
		
	# If method is equal
	elif [ $METHOD = 1 ]
	then
		echo "Equal mixing method selected"
		cd $OUTPUT_DIR/re-mix
		COMBINED_LINES=$(wc -l < $TREATMENT_NAME_BASE.bootstrap_$BOOTSTRAP_COUNT.tmp)
		FIRST_LINES=$(( COMBINED_LINES / 2 ))
		LAST_LINES=$(( COMBINED_LINES - FIRST_LINES ))
	fi
	
	echo "Re-mixing..."
	shuf -o $TREATMENT_NAME_BASE.bootstrap_$BOOTSTRAP_COUNT.tmp < $TREATMENT_NAME_BASE.bootstrap_$BOOTSTRAP_COUNT.tmp
	echo "Check!"
	
	echo "Creating re-mixed treatment library"
	head -n $FIRST_LINES $TREATMENT_NAME_BASE.bootstrap_$BOOTSTRAP_COUNT.tmp > $TREATMENT_NAME_BASE.bootstrap_$BOOTSTRAP_COUNT.bed
	echo "Check!"
	
	echo "Creating re-mixed control library"
	tail -n $LAST_LINES  $TREATMENT_NAME_BASE.bootstrap_$BOOTSTRAP_COUNT.tmp > $CONTROL_NAME_BASE.bootstrap_$BOOTSTRAP_COUNT.bed
	echo "Check!"
	
	echo "Deleting temporary files"
	rm $TREATMENT_NAME_BASE.bootstrap_$BOOTSTRAP_COUNT.tmp
	echo "Check!"
	echo "Completed re-mix #$BOOTSTRAP_COUNT"
done

cd $INPUT_DIR

end_time=`date +%s`
echo RECAP RE-MIX execution time: `expr $end_time - $start_time`s.

#################################################################
######################## End Script Here ########################
}


#################### Begin Options and Usage ####################
# Print usage
usage() {
  echo -n "
  [Input directory]  [Treatment file] [Control file] 
  [Output directory] [Re-mix method]  [Bootstrap]

 ${bold}Options:${normal}
  -i, --input 	    Input file directory
  -t, --treatment   Treatment file
  -c, --control     Control file
  -o, --output      Output file directory
  -m, --method      Method of re-mixing (equal) or (unequal)
  -b, --bootstrap   Number of re-mixes
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
    	-t|--treatment) shift; TREATMENT_NAME=${1} ;;
    	-c|--control)   shift; CONTROL_NAME=${1} ;;
    	-o|--output)    shift; OUTPUT_DIR=${1} ;;
    	-m|--method)    shift 
    					if [ ${1} = 'equal' ]
    					then
							METHOD=1
							METHOD_NAME="Equal library sizes"
						elif [ ${1} = 'unequal' ]
						then
							METHOD=0
							METHOD_NAME="Default library sizes"
						else
							echo "ERROR: Invalid option"
							exit 1
						fi ;;
		-b|--bootstrap) shift; BOOTSTRAP=${1} ;;
		-h|--help)      usage >&2; exit 0 ;;
		*)              echo "ERROR: Bad argument ${1}" ; exit 1 ;;
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
mainRemix
# ===============================================================
