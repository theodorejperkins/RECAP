#!/usr/bin/perl
# ==================================================================
# RECAP Re-Mix
# RECAP is a wrapper algorithm that resamples ChIP-seq and control
# data to estimate and control for biases built into peak calling
# algorithms.
# The purpose of this script is to recalibrate p-values based on
# the original ChIP-seq/control and re-mixed ChIP-seq/control
# data. A detailed explanation of the algorithm can be found
# under PUBLICATIONS.
#
# HISTORY:
#   02/01/2018 - v1.0.0 - First Creation
#   02/08/2018 - v1.0.1 - Outputs header and new column for RECAP p-values
#   02/12/2018 - v1.0.2 - Binary search instead of first index search
#   03/11/2018 - v1.0.3 - Filters off downregulated diffReps p-values
#   03/14/2018 - v1.0.4 - Calculates BH values for RECAP r-values
#   03/22/2018 - v1.0.5 - Calculates the LFDR binned by order of magitude
#   05/12/2018 - v1.0.6 - LFDR binned by half decade
#   05/28/2018 - v1.0.7 - Calculates bootstrapped RECAP p-values
#   05/29/2018 - v1.0.8 - Saves output to separate bootstrapped files
#   08/24/2018 - v1.1.0 - Major bugfix to RECAP calculation with tied p-values
#   08/28/2018 - v1.1.1 - Relative path fix when using ./
#   11/18/2018 - v1.2.0 - Linear interpolation of RECAP p-values to remove duplicates
#   11/20/2018 - v1.2.1 - Linear interpolation changes
#   11/21/2018 - v1.2.2 - More linear interpolation changes
#   12/24/2018 - v1.2.3 - Funky distribution tail fixer (pseudocount re-mixed p-value)
#   01/14/2019 - v1.2.4 - Remove all LFDR and no bootstrap folder
#   01/16/2019 - v1.2.5 - Better instructions
#
# CREDITS:
# RECAP was developed by Justin G. Chitpin, Aseel Awdeh, and Theodore J. Perkins.
# Development of RECAP was carried out at the Ottawa Hospital
# Research Institute in the Perkins Lab.
#
# PUBLICATIONS:
# If you use RECAP, please cite the following paper:
# <INSERT PUBLICATION HERE>
#
# QUESTIONS:
# Please contact tperkins@ohri.ca
# ==================================================================


# ==================================================================
# Version
our $version = "1.2.5";

# Modules
use warnings;
use strict;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use List::BinarySearch qw(binsearch_pos);
use List::Util qw(min max);
use List::MoreUtils qw(uniq);
use Math::Utils qw(ceil floor);
use autodie;
$| = 1;
# ==================================================================

######################## Begin Script Here #########################
####################################################################

# Start job time
BEGIN { our $start_run = time(); }

# Variables
my @remixList;    # Holds the names of bootstrapped re-mixed files
my @fileOriginal; # Holds contents of original file
my @fileRemix;    # Holds contents of re-mixed file
my @fileRemix2;   # Holds contents of multiple re-mixed files
my @columnsOriginal; # Holds contents of original file split by delimiter
my @columnsRemix;    # Holds contents of original file split by delimiter
my @pvalOrig;    # Original p-values
my @pvalRemix;   # Re-mixed p-values
my @RECAP;       # Re-calibrated r-values
my @fileHeader;  # File header for output file
my $delimOutput; # Delimiter type passed
my @indexFDR;    # Index of FDR values corresponding to original file
my @joinFDR;     # FDR values indexed by $indexFDR

# Arguments
GetOptions(
  'dirOrig=s'    => \my $dirOriginal,
  'nameOrig=s'   => \my $nameOriginal,
  'dirRemix=s'   => \my $dirRemix,
  'nameRemix=s'  => \my $nameRemix,
  'dirOutput=s'  => \my $dirOutput,
  'nameOutput=s' => \my $nameOutput,
  'bootstrap=i'  => \my $bootstrap,
  'header=i'     => \my $header,
  'pvalCol=i'    => \my $pvalCol,
  'delim=s'      => \my $delim,
  'software=s'   => \my $software,
  'help'         => \my $help,
) or die "ERROR: Invalid options passed to $0\n"; 
# List parameters and check if all arguments filled
if ( $help || ! defined $dirOriginal
           || ! defined $nameOriginal
           || ! defined $dirRemix
           || ! defined $dirOutput
           || ! defined $nameOutput
           || ! defined $nameRemix
           || ! defined $bootstrap
           || ! defined $header
           || ! defined $pvalCol
           || ! defined $delim
           || ! defined $software) {
  print BOLD, "  USAGE:\n", RESET;
  print "  --dirOrig\tOriginal file directory (absolute)\n";
  print "  --nameOrig\tOriginal file name\n";
  print "  --dirRemix\tRe-mixed file directory (absolute)\n";
  print "  --nameRemix\tRe-mixed file name without '.bootstrap_#.bed\n";
  print "  --dirOutput\tOutput file directory (absolute)\n";
  print "  --nameOutput\tOutput file name\n";
  print "  --bootstrap\tNumber of re-mixing procedures\n";
  print "  --header\tHeader lines in files\n";
  print "  --pvalCol\tColumn containing p-values\n";
  print "  --delim\tDelimiter type of file (t)ab or (c)omma\n";
  print "  --software\tPeak caller (M)ACS2 or (D)iffReps or (O)ther\n\n";
  print BOLD, "  OPTIONS:\n", RESET;
  print "  --help\tDisplay this help and exit\n";
  exit;
}

# Standard output
print "##################################################\n";
print "############       RECAP v$version       ############\n";
print "##################################################\n";
print "\n";
print "Input directory  (Original; absolute path): $dirOriginal\n";
print "Input file name  (Original): $nameOriginal\n";
print "Input directory  (Re-mixed; absolute path): $dirRemix\n";
print "Input file name  (Re-mixed): $nameRemix\n";
print "Output directory (Absolute path): $dirOutput\n";
print "Output file name: $nameOutput\n";
print "Bootstrap number: $bootstrap\n";
print "Header number   : $header\n";
print "p-value column  : $pvalCol\n";
print "Delimiter type  : $delim\n";
print "Peak caller type: $software\n";

# Input validation
if ( ! -d $dirOriginal ) {
	die "ERROR: Directory of original file does not exist\n";
}

chdir $dirOriginal;
if ( ! -e $nameOriginal ) {
	die "ERROR: Original file does not exist\n";
}

if ( ! -d $dirRemix ) {
	die "ERROR: Directory of re-mixed file does not exist\n";
}

{
	my $numFiles;
	opendir (my $DIR, $dirRemix ); 
	# goatse operator
	$numFiles = () = grep { /$nameRemix/ && /bootstrap/ } readdir $DIR;
	closedir $DIR; 
	if ( $numFiles == 0 ) {
		die "ERROR: Re-mixed file does not exist\n";
	}
	if ( $bootstrap > $numFiles ) {
		die "ERROR: Bootstrap number exceeds number of re-mixed files\n";
	} 
}

if ( ! -d $dirOutput ) {
	die "ERROR: Directory for output does not exist\n";
}

if ( ! $bootstrap == [0-9]*([0-9]) ) {
	die "ERROR: Specify a whole number for the bootstrap\n";
}

if ( ! $header == [0-9]*([0-9]) ) {
	die "ERROR: Specify a whole number for the header\n";
}

if ( ! $pvalCol == [0-9]*([0-9]) ) {
	die "ERROR: Specify a whole number for the p-value column\n";
} else {
	$pvalCol--; # p-value column must be n-1 for split() to work
}

$delim = lc $delim;
if ( $delim =~ /[^ct]/ ) {
	die "ERROR: Specify 'c' for comma or 't' for tab\n";
} else {
	if ( $delim eq "c" ) {
		$delim = qr/,/;
		$delimOutput = ",";
	} elsif ( $delim eq "t" ) {
		$delim = qr/\t/;
		$delimOutput = "\t";
	}
}

$software = uc $software;
if ( $software =~ /[^MDO]/ ) {
	die "ERROR: Specify 'M' for MACS, 'D' for diffReps, 'O' for other\n";
}

# Reading original summary file
print "\n";
print "Reading original summary file\n";
chdir($dirOriginal);
open(my $importOrig, '<', $nameOriginal);
while(<$importOrig>) {
	next unless $. > $header;
	push @fileOriginal, $_;
}
close $importOrig;
print "Check!\n";

# Reading re-mixed summary file(s)
print "Reading multiple re-mixed summary files\n";
chdir($dirRemix);

# Scanning names of all re-mixed summary files in ascending order 
opendir(my $DIR, $dirRemix);
@remixList = grep { /$nameRemix/ && /bootstrap/ } readdir $DIR;
closedir $DIR;
@remixList = sort { $a cmp $b } @remixList;

# Pushing all re-mixed datasets into one array
for ( my $i=0; $i<$bootstrap; $i++ ) {
	open(my $importRemix, '<', $remixList[$i]);
	while(<$importRemix>) {
		next unless $. > $header;
		push @fileRemix2, $_;
	}
	close $importRemix;
	push @fileRemix, @fileRemix2;
}
print "Check!\n";

# Removing downregulated p-values if diffReps summary file
if ( $software eq "D" ) {
	print "Filtering off downregulated diffReps p-values from original/re-mixed files\n";
	@fileOriginal = grep( /Up/, @fileOriginal );
	@fileRemix = grep( /Up/, @fileRemix );
}

# Extracting p-value columns from summary files
# Exponentiate p-values if MACS summary files
print "Scanning p-value column of original file\n";
foreach my $line (@fileOriginal) {
	@columnsOriginal = split(/$delim/, "$line");
	if ( $software eq "M" ) {
		push @pvalOrig, 10**-($columnsOriginal[$pvalCol]);
	} else {
		push @pvalOrig, $columnsOriginal[$pvalCol];
	}
}
print "Check!\n";
print "Scanning p-value column of re-mixed file\n";
foreach my $line (@fileRemix) {
	@columnsRemix = split(/$delim/, "$line");
	if ( $software eq "M" ) {
		push @pvalRemix, 10**-($columnsRemix[$pvalCol]);
	} else {
		push @pvalRemix, $columnsRemix[$pvalCol];
	}
}
print "Check!\n";

## RECAP algorithm
# Implementing RECAP algorithm with binary search
{
	# Sort re-mixed p-values
	print "Sorting re-mixed p-values\n";
	my @sorted_pvalRemix = @pvalRemix;
	push @sorted_pvalRemix, 0;
	@sorted_pvalRemix = sort { $b <=> $a } @sorted_pvalRemix;
	print "Check!\n";

	print "Recalibrating via RECAP procedure\n";
	for my $idx ( 0 .. $#pvalOrig ) {
		# Number of re-mixed p-values equal or less than original p
		my $matching = binsearch_pos { $b <=> $a } $pvalOrig[$idx], @sorted_pvalRemix;
		$RECAP[$idx] = ((scalar @sorted_pvalRemix) - $matching) / scalar @sorted_pvalRemix;
	}
	print "Check!\n";
}

# Linear interpolation of RECAP p-values
print "Linear interpolation of RECAP p-values based on original and re-mixed p-values\n";
{
	# Sorting re-mixed p-values and saving a vector of unique values
	my @sorted_pvalRemix = @pvalRemix;
	push @sorted_pvalRemix, 0;
	@sorted_pvalRemix = sort { $b <=> $a } @sorted_pvalRemix;
	my @sorted_pvalRemix2 = @pvalRemix;
	push @sorted_pvalRemix2, 0;
	@sorted_pvalRemix2 = sort { $a <=> $b } @sorted_pvalRemix;

	#my @sorted_pvalRemix = sort { $b <=> $a } @pvalRemix;
	#my @sorted_pvalRemix2 = sort { $a <=> $b } @pvalRemix;
	my @sorted_pvalRemix_uniq = uniq @sorted_pvalRemix2;

	# Initialize linear interpolated RECAP p-values
	my @liRECAP;

	# Linear interpolation code begins here
	for my $idx ( 0 .. $#pvalOrig ) {
		my $matching = binsearch_pos { $a <=> $b } $pvalOrig[$idx], @sorted_pvalRemix_uniq;
		my $r2 = $sorted_pvalRemix_uniq[$matching];   # the bigger re-mixed p-value vs $pvalOrig[$idx]
		my $r1 = $sorted_pvalRemix_uniq[$matching-1]; # the smaller re-mixed p-value vs $pvalOrig[$idx]
		if ( $pvalOrig[$idx] < $sorted_pvalRemix_uniq[0] ) {
			my $r1Matching = binsearch_pos { $b <=> $a } $sorted_pvalRemix_uniq[0], @sorted_pvalRemix;
			my $r1RECAP = (( scalar @sorted_pvalRemix) - $r1Matching) / scalar @sorted_pvalRemix;
			my $temp = $pvalOrig[$idx]/$sorted_pvalRemix_uniq[0] * $r1RECAP;
			push @liRECAP, $temp;
		} elsif ( defined $r2 && $pvalOrig[$idx] == $r2) {
			push @liRECAP, $RECAP[$idx];
		} elsif ( defined $r2 && $r1 <= $pvalOrig[$idx] && $pvalOrig[$idx] <= $r2 ) {
			my $r2Matching = binsearch_pos { $b <=> $a } $r2, @sorted_pvalRemix;
			my $r2RECAP = (( scalar @sorted_pvalRemix) - $r2Matching) / scalar @sorted_pvalRemix;
			my $r1Matching = binsearch_pos { $b <=> $a } $r1, @sorted_pvalRemix;
			my $r1RECAP = (( scalar @sorted_pvalRemix) - $r1Matching) / scalar @sorted_pvalRemix;
			my $temp = ($pvalOrig[$idx] - $r1)/($r2-$r1) * ($r2RECAP - $r1RECAP) + $r1RECAP;
			push @liRECAP, $temp;
		} elsif ( $pvalOrig[$idx] > $sorted_pvalRemix_uniq[$#sorted_pvalRemix_uniq] ) {
			push @liRECAP, 1;
		} else {
			print("$pvalOrig[$idx] > $sorted_pvalRemix_uniq[$#sorted_pvalRemix_uniq] || $r1 <= $pvalOrig[$idx] <= $r2");
			die("Linear interpolation failed. Check line 313 in the source code.");
		}
	}
	# Overwrite original RECAP p-values with the interpolated ones
	@RECAP = @liRECAP;
	print "Check!\n";
}

## Benjamini-Hochberg adjustment of RECAP r-values
print "FDR-adjusting RECAP r-values\n";

# Join RECAP values with its original index
@indexFDR = map{ $_ } 1 .. scalar @RECAP;
@joinFDR = map[ $indexFDR[$_], $RECAP[$_] ], 0 .. $#indexFDR;

# Sort 2D array by RECAP r-value
@joinFDR = sort{ $a->[1] <=> $b->[1] } @joinFDR;

# Implementing BH, allowing r-values of the same rank
{
	my $rank=0;
	for ( my $i=0; $i<scalar @joinFDR; $i++ ) {
		unless ( $i == scalar @joinFDR && $joinFDR[$i] == $joinFDR[$i+1] ) {
			$rank++;
		}
		$joinFDR[$i]->[1] = $joinFDR[$rank-1][1] / ($rank / scalar @joinFDR);
	}
}

# Adjusting FDR for the next higher raw r-value if it is smaller than r-value times m/i
for ( my $i=$#joinFDR; $i>0; $i-- ) {
	if ( $joinFDR[$i][1] < $joinFDR[$i-1][1] ) {
		$joinFDR[$i-1]->[1] = $joinFDR[$i][1];
	}
}

# Re-sort 2D array of FDR r-values by original index
@joinFDR = sort { $a->[0] <=> $b->[0] } @joinFDR;
# Only retain column of FDR r-values
@joinFDR = map{ $$_[1] } @joinFDR;

# Save recalibrated p-values to file
print "Saving RECAP summary file to output\n";

# Saving header from original file and adding RECAP BH(RECAP) column headers
if ( $header != 0 ) {
  chdir($dirOriginal);
  open(my $importOrig, '<', $nameOriginal);
  while(<$importOrig>) {
    next unless $. <= $header;
    push @fileHeader, $_;
  }
chomp($fileHeader[$header-1]);
push @fileHeader, "${delimOutput}RECAP${delimOutput}BH(RECAP)${delimOutput}\n";
}

# Create a new folder called "bootstrap_#" and save output to that folder
chdir($dirOutput);
open ( my $exportRECAP, ">", $nameOutput );
print $exportRECAP join "", @fileHeader;

# Add column of r-values to input, original summary file
# Add column of Benjamini-Hochberg corrected r-value
for my $idx ( 0 .. $#fileOriginal ) {
  @columnsOriginal = split(/$delim/, "$fileOriginal[$idx]");
  chomp(@columnsOriginal);
  push @columnsOriginal, $RECAP[$idx], $joinFDR[$idx], "\n";
  print $exportRECAP join $delimOutput, @columnsOriginal;
}
close $exportRECAP;
print "Check!\n";

# Job time
my $end_run = time();
my $run_time = $end_run - our $start_run;
print "Job took $run_time seconds\n";
