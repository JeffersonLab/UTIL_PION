#!/bin/bash

# 10/08/21 - Stephen JD Kay - University of Regina
# This script sets up the relevant sym links for running analysis scripts
# This version is specific for cdaq

ANALYSISPATH="/home/cdaq/pionLT-2021"
OUTPUTPATH="/net/cdaq/cdaql1data/cdaq/hallc-online-pionLT"

# Make the sym links
echo "Sym links will now be made for - "
echo ""
echo "${OUTPUTPATH}/ROOTfiles/ - ROOTfiles"
echo "${OUTPUTPATH}/OUTPUT/ - OUTPUT"
echo "${OUTPUTPATH}/REPORT_OUTPUT/ - REPORT_OUTPUT"
echo "${OUTPUTPATH}/HISTOGRAMS/ - HISTOGRAMS"
echo""
read -p "Proceed with sym link setup setup? <Y/N> " prompt2
if [[ $prompt2 == "y" || $prompt2 == "Y" || $prompt2 == "yes" || $prompt2 == "Yes" ]]; then
    echo "Creating sym links"
else 
    echo "Please update script with desired changes to sym links and re-run if desired"
    exit 1
fi

# Each loop checks if the link exists, if it doesn't, make it
# If it DOES, check it's not broken, if broken, replace, if not, just print that it exists
if [ ! -L "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/ROOTfiles" ]; then
    ln -s "${OUTPUTPATH}/ROOTfiles/" "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/ROOTfiles"
elif [ -L "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/ROOTfiles" ]; then
    if [ ! -e "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/ROOTfiles" ]; then
	echo "ROOTfiles sym link exits but is broken, replacing"
	rm "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/ROOTfiles"
	ln -s "${OUTPUTPATH}/ROOTfiles/" "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/ROOTfiles"
    else echo "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/ROOTfiles sym link already exists and not broken"
    fi
fi

if [ ! -L "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/OUTPUT" ]; then
    ln -s "${OUTPUTPATH}/OUTPUT/" "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/OUTPUT"
elif [ -L "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/OUTPUT" ]; then
    if [ ! -e "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/OUTPUT" ]; then
	echo "OUTPUT sym link exits but is broken, replacing"
	rm "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/OUTPUT"
	ln -s "${OUTPUTPATH}/OUTPUT/" "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/OUTPUT"
    else echo "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/OUTPUT sym link already exists and not broken"
    fi
fi

if [ ! -L "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/REPORT_OUTPUT" ]; then
    ln -s "${OUTPUTPATH}/REPORT_OUTPUT/" "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/REPORT_OUTPUT"
elif [ -L "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/REPORT_OUTPUT" ]; then
    if [ ! -e "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/REPORT_OUTPUT" ]; then
	echo "REPORT_OUTPUT sym link exits but is broken, replacing"
	rm "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/REPORT_OUTPUT"
	ln -s "${OUTPUTPATH}/REPORT_OUTPUT/" "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/REPORT_OUTPUT"
    else echo "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/REPORT_OUTPUT sym link already exists and not broken"
    fi
fi

if [ ! -L "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/HISTOGRAMS" ]; then
    ln -s "${OUTPUTPATH}/HISTOGRAMS/" "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/HISTOGRAMS"
elif [ -L "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/HISTOGRAMS" ]; then
    if [ ! -e "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/HISTOGRAMS" ]; then
	echo "HISTOGRAMS sym link exits but is broken, replacing"
	rm "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/HISTOGRAMS"
	ln -s "${OUTPUTPATH}/HISTOGRAMS/" "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/HISTOGRAMS"
    else echo "${ANALYSISPATH}/hallc_replay_lt/UTIL_PION/HISTOGRAMS sym link already exists and not broken"
    fi
fi

echo ""
echo "Directories and sym links for UTIL_PION now setup"

exit 0
