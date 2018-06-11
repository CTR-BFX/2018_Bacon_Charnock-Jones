#!/bin/bash

#------------------------------------------------------------------------------
# Demultiplex Raw Fastq files for DropSeq 
# Wrapper script for BBTools demuxbyname.sh 
# Available from https://sourceforge.net/projects/bbmap/
#
# To accompany:
# Bacon et al (2018) Coup of T: Finding missing immune cell subtypes in growth restricted neonates 
# using novel single-cell sequencing analysis
#
# Link to publication
# TO ADD ONCE AVAILABLE
#
# Script available from:
# https://github.com/CTR-BFX/2018_Bacon_Charnock-Jones
#
# Analysis Performed by Russell S. Hamilton
# CTR Bioinformatics 
# Centre for Trophoblast Reseach, University of Cambridge, UK
# Copyright Russell S. Hamilton (rsh46@cam.ac.uk)
#
#------------------------------------------------------------------------------
# License Information
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------

set -e

usage() { echo "Usage: $0 -e <SLX-10282>" 1>&2; exit 1; }


NEXTERA=("N701:TAAGGCGA"
         "N702:CGTACTAG"
         "N703:AGGCAGAA"
         "N704:TCCTGAGC"
         "N705:GGACTCCT"
         "N706:TAGGCATG"
         "N707:CTCTCTAC"
         "N710:CGAGGCTG"
         "N711:AAGAGGCA"
         "N712:GTAGAGGA"
         "N714:GCTCATGA" 
         "N715:ATCTCAGG"
         "N716:ACTCGCTA"
         "N718:GGAGCTAC"
         "N719:GCGTAGTA"
         "N720:CGGAGCCT"
         "N721:TACGCTGC"
         "N722:ATGCGCAG"
         "N723:TAGCGCTC"
         "N724:ACTGAGCG"
         "N726:CCTAAGAC"
         "N727:CGATCAGT"
         "N728:TGCAGCTA"
         "N729:TCGACGTC"
        )

while getopts ":e:i:n:" o; do
    case "${o}" in
        e)
            e=${OPTARG}
            ;; 
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${e}" ]; then
    usage
fi

IDXs=""

for nextera in "${NEXTERA[@]}" ; 
  do
     KEY="${nextera%%:*}"
     VALUE="${nextera##*:}"
     IDXs="$IDXs,$VALUE"

     #echo $nextera ":::" $KEY ":::" $VALUE
  done

IDXstring=`echo $IDXs | sed 's/^,//g'`


echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "+ Running: $0"
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "+ Contact: Russell S. Hamilton, 2016, rsh46@cam.ac.uk"
echo "+ Experiment = " ${e}
echo "+ NEXTERA    = " ${IDXstring}

# Has been NoIndex and UnspecifiedIndex
#IDXSTR="NoIndex"
IDXSTR="UnspecifiedIndex"

FILE=${e}.${IDXSTR}.*.s_[1-8].r_1.*.gz
FILESTEM=`echo ${FILE} | sed 's/.s_[1-8].r_1.*//g'`
LANE=`ls -1 ${FILE} | sed 's/.r_.*.fq.gz//' | sed 's/.*s_//' `



declare -i TOTAL1
declare -i TOTAL2
TOTAL1=`gunzip -c ${e}.${IDXSTR}.*.s_${LANE}.r_1.*.gz | grep -c "^@[A-Z]"`
TOTAL2=`gunzip -c ${e}.${IDXSTR}.*.s_${LANE}.r_2.*.gz | grep -c "^@[A-Z]"`

# Test only for speed up as the above step is very slow
#TOTAL1=228843716
#TOTAL2=228843716


echo "+ Reads in R1 fq.gz = " $TOTAL1
echo "+ Reads in R2 fq.gz = " $TOTAL2


demuxbyname.sh in=${FILESTEM}.s_${LANE}.r_#.fq.gz out=${e}.%.s_${LANE}.r_#.fq.gz prefixmode=f substringmode=t names=${IDXstring}

for nextera in "${NEXTERA[@]}" ; 
  do
     KEY="${nextera%%:*}"
     VALUE="${nextera##*:}"
     IDXs="$IDXs,$VALUE"

     R1FileOld=`echo ${e}"."$VALUE"."s_${LANE}.r_1.fq.gz`
     R2FileOld=`echo ${e}"."$VALUE"."s_${LANE}.r_2.fq.gz`

     R1FileNew=`echo ${e}"."$VALUE"."$KEY"."s_${LANE}.r_1.fq.gz`
     R2FileNew=`echo ${e}"."$VALUE"."$KEY"."s_${LANE}.r_2.fq.gz`

     mv $R1FileOld $R1FileNew
     mv $R2FileOld $R2FileNew

     declare -i CNT1
     declare -i CNT2
     CNT1=`gunzip -c $R1FileNew | grep -c "^@[A-Z][A-Z0-9]*:[0-9]*:"`
     CNT2=`gunzip -c $R2FileNew | grep -c "^@[A-Z][A-Z0-9]*:[0-9]*"`

     echo "+ CNT1 = " ${CNT1}
     echo "+ CNT2 = " ${CNT2}

     PERC1=`echo "scale=4; 100*($CNT1/$TOTAL1)" | bc`
     PERC2=`echo "scale=4; 100*($CNT2/$TOTAL2)" | bc`

     echo "+" $R1FileNew "=" ${PERC1} "=" ${CNT1}
     echo "+" $R2FileNew "=" ${PERC2} "=" ${CNT2}

     minperc=0.99
     if [[ "${PERC1}" -le "$minperc" ]]; then
     	echo "Remove $R1FileNew?"
     fi
     
     if [[ "${PERC2}" -le "$minperc" ]]; then
     	echo "Remove $R2FileNew?"
     fi 	
  done

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "+ Finished: $0"
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

#
# FIN
#
