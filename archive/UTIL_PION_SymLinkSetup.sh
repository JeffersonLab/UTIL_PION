#!/bin/bash

# 24/02/21 - Stephen JD Kay - University of Regina
# This script sets up the relevant sim links for the UTIL_PION scripts
# The directories should already exists if you ran the script in hallc_replay_lt
# It makes the relevant directories in /volatile if they do not exist

echo "Beginning setup of folders and symlinks"
# Set path depending upon hostname. Change or add more as needed  
if [[ "${HOSTNAME}" = *"farm"* ]]; then 
    GROUPPATH="/group/c-pionlt"
    USERPATH="${GROUPPATH}/USERS/${USER}"
    VOLATILEPATH="/volatile/hallc/c-pionlt"
    cd "${USERPATH}/hallc_replay_lt"
else echo "Host not recognised, please add relevant pathing for hostname to the script and re-run"
fi

# Next, need to make a load of directories, check they exist, if not, make em!
if [ ! -d "${VOLATILEPATH}/${USER}" ]; then
    mkdir "${VOLATILEPATH}/${USER}"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/ROOTfiles" ]; then
    mkdir "${VOLATILEPATH}/${USER}/ROOTfiles"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/OUTPUT" ]; then
    mkdir "${VOLATILEPATH}/${USER}/OUTPUT"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/REPORT_OUTPUT" ]; then
    mkdir "${VOLATILEPATH}/${USER}/REPORT_OUTPUT"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/HISTOGRAM" ]; then
    mkdir "${VOLATILEPATH}/${USER}/HISTOGRAMS"
fi

# Now, check if subdirectories exist for ROOTfiles, if not, make them
# Analysis subdirectories
if [ ! -d "${VOLATILEPATH}/${USER}/ROOTfiles/Analysis" ]; then
    mkdir "${VOLATILEPATH}/${USER}/ROOTfiles/Analysis"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/ROOTfiles/Analysis/HeeP" ]; then
    mkdir "${VOLATILEPATH}/${USER}/ROOTfiles/Analysis/HeeP"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/ROOTfiles/Analysis/Lumi" ]; then
    mkdir "${VOLATILEPATH}/${USER}/ROOTfiles/Analysis/Lumi"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/ROOTfiles/Analysis/PionLT" ]; then
    mkdir "${VOLATILEPATH}/${USER}/ROOTfiles/Analysis/PionLT"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/ROOTfiles/Analysis/PID" ]; then
    mkdir "${VOLATILEPATH}/${USER}/ROOTfiles/Analysis/PID"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/ROOTfiles/Analysis/Optics" ]; then
    mkdir "${VOLATILEPATH}/${USER}/ROOTfiles/Analysis/Optics"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/ROOTfiles/Analysis/General" ]; then
    mkdir "${VOLATILEPATH}/${USER}/ROOTfiles/Analysis/General"
fi
# Scalers ROOTfiles directory
if [ ! -d "${VOLATILEPATH}/${USER}/ROOTfiles/Scalers" ]; then
    mkdir "${VOLATILEPATH}/${USER}/ROOTfiles/Scalers"
fi
# Now, check if subdirectories exist for OUTPUT, if not, make them
# Analysis subdirectories
if [ ! -d "${VOLATILEPATH}/${USER}/OUTPUT/Analysis" ]; then
    mkdir "${VOLATILEPATH}/${USER}/OUTPUT/Analysis"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/OUTPUT/Analysis/HeeP" ]; then
    mkdir "${VOLATILEPATH}/${USER}/OUTPUT/Analysis/HeeP"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/OUTPUT/Analysis/Lumi" ]; then
    mkdir "${VOLATILEPATH}/${USER}/OUTPUT/Analysis/Lumi"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/OUTPUT/Analysis/PionLT" ]; then
    mkdir "${VOLATILEPATH}/${USER}/OUTPUT/Analysis/PionLT"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/OUTPUT/Analysis/PID" ]; then
    mkdir "${VOLATILEPATH}/${USER}/OUTPUT/Analysis/PID"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/OUTPUT/Analysis/Optics" ]; then
    mkdir "${VOLATILEPATH}/${USER}/OUTPUT/Analysis/Optics"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/OUTPUT/Analysis/General" ]; then
    mkdir "${VOLATILEPATH}/${USER}/OUTPUT/Analysis/General"
fi
# Scalers OUTPUT directory
if [ ! -d "${VOLATILEPATH}/${USER}/OUTPUT/Scalers" ]; then
    mkdir "${VOLATILEPATH}/${USER}/OUTPUT/Scalers"
fi
# Now, check if subdirectories exist for REPORT_OUTPUT, if not, make them
# Analysis subdirectories
if [ ! -d "${VOLATILEPATH}/${USER}/REPORT_OUTPUT/Analysis" ]; then
    mkdir "${VOLATILEPATH}/${USER}/REPORT_OUTPUT/Analysis"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/REPORT_OUTPUT/Analysis/HeeP" ]; then
    mkdir "${VOLATILEPATH}/${USER}/REPORT_OUTPUT/Analysis/HeeP"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/REPORT_OUTPUT/Analysis/Lumi" ]; then
    mkdir "${VOLATILEPATH}/${USER}/REPORT_OUTPUT/Analysis/Lumi"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/REPORT_OUTPUT/Analysis/PionLT" ]; then
    mkdir "${VOLATILEPATH}/${USER}/REPORT_OUTPUT/Analysis/PionLT"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/REPORT_OUTPUT/Analysis/PID" ]; then
    mkdir "${VOLATILEPATH}/${USER}/REPORT_OUTPUT/Analysis/PID"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/REPORT_OUTPUT/Analysis/Optics" ]; then
    mkdir "${VOLATILEPATH}/${USER}/REPORT_OUTPUT/Analysis/Optics"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/REPORT_OUTPUT/Analysis/General" ]; then
    mkdir "${VOLATILEPATH}/${USER}/REPORT_OUTPUT/Analysis/General"
fi
# Scalers REPORT_OUTPUT directory
if [ ! -d "${VOLATILEPATH}/${USER}/REPORT_OUTPUT/Scalers" ]; then
    mkdir "${VOLATILEPATH}/${USER}/REPORT_OUTPUT/Scalers"
fi
# Now, check if subdirectories exist for HISTOGRAMS, if not, make them
# Analysis subdirectories
if [ ! -d "${VOLATILEPATH}/${USER}/HISTOGRAMS/Analysis" ]; then
    mkdir "${VOLATILEPATH}/${USER}/HISTOGRAMS/Analysis"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/HISTOGRAMS/Analysis/HeeP" ]; then
    mkdir "${VOLATILEPATH}/${USER}/HISTOGRAMS/Analysis/HeeP"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/HISTOGRAMS/Analysis/Lumi" ]; then
    mkdir "${VOLATILEPATH}/${USER}/HISTOGRAMS/Analysis/Lumi"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/HISTOGRAMS/Analysis/PionLT" ]; then
    mkdir "${VOLATILEPATH}/${USER}/HISTOGRAMS/Analysis/PionLT"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/HISTOGRAMS/Analysis/PID" ]; then
    mkdir "${VOLATILEPATH}/${USER}/HISTOGRAMS/Analysis/PID"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/HISTOGRAMS/Analysis/Optics" ]; then
    mkdir "${VOLATILEPATH}/${USER}/HISTOGRAMS/Analysis/Optics"
fi
if [ ! -d "${VOLATILEPATH}/${USER}/HISTOGRAMS/Analysis/General" ]; then
    mkdir "${VOLATILEPATH}/${USER}/HISTOGRAMS/Analysis/General"
fi
# Scalers HISTOGRAMS directory
if [ ! -d "${VOLATILEPATH}/${USER}/HISTOGRAMS/Scalers" ]; then
    mkdir "${VOLATILEPATH}/${USER}/HISTOGRAMS/Scalers"
fi
# Directories now created, make the sym links
echo "Directories have been created in volatile if they did not already exists, sym links will now be made for - "
echo ""
echo "${VOLATILEPATH}/${USER}/ROOTfiles/ - ROOTfiles"
echo "${VOLATILEPATH}/${USER}/OUTPUT/ - OUTPUT"
echo "${VOLATILEPATH}/${USER}/REPORT_OUTPUT/ - REPORT_OUTPUT"
echo "${VOLATILEPATH}/${USER}/HISTOGRAMS/ - HISTOGRAMS"
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

if [ ! -L "${USERPATH}/hallc_replay_lt/UTIL_PION/ROOTfiles" ]; then
    ln -s "${VOLATILEPATH}/${USER}/ROOTfiles/" "${USERPATH}/hallc_replay_lt/UTIL_PION/ROOTfiles"
elif [ -L "${USERPATH}/hallc_replay_lt/UTIL_PION/ROOTfiles" ]; then
    if [ ! -e "${USERPATH}/hallc_replay_lt/UTIL_PION/ROOTfiles" ]; then
	echo "ROOTfiles sym link exits but is broken, replacing"
	rm "${USERPATH}/hallc_replay_lt/UTIL_PION/ROOTfiles"
	ln -s "${VOLATILEPATH}/${USER}/ROOTfiles/" "${USERPATH}/hallc_replay_lt/UTIL_PION/ROOTfiles"
    else echo "${USERPATH}/hallc_replay_lt/UTIL_PION/ROOTfiles sym link already exists and not broken"
    fi
fi

if [ ! -L "${USERPATH}/hallc_replay_lt/UTIL_PION/OUTPUT" ]; then
    ln -s "${VOLATILEPATH}/${USER}/OUTPUT/" "${USERPATH}/hallc_replay_lt/UTIL_PION/OUTPUT"
elif [ -L "${USERPATH}/hallc_replay_lt/UTIL_PION/OUTPUT" ]; then
    if [ ! -e "${USERPATH}/hallc_replay_lt/UTIL_PION/OUTPUT" ]; then
	echo "OUTPUT sym link exits but is broken, replacing"
	rm "${USERPATH}/hallc_replay_lt/UTIL_PION/OUTPUT"
	ln -s "${VOLATILEPATH}/${USER}/OUTPUT/" "${USERPATH}/hallc_replay_lt/UTIL_PION/OUTPUT"
    else echo "${USERPATH}/hallc_replay_lt/UTIL_PION/OUTPUT sym link already exists and not broken"
    fi
fi

if [ ! -L "${USERPATH}/hallc_replay_lt/UTIL_PION/REPORT_OUTPUT" ]; then
    ln -s "${VOLATILEPATH}/${USER}/REPORT_OUTPUT/" "${USERPATH}/hallc_replay_lt/UTIL_PION/REPORT_OUTPUT"
elif [ -L "${USERPATH}/hallc_replay_lt/UTIL_PION/REPORT_OUTPUT" ]; then
    if [ ! -e "${USERPATH}/hallc_replay_lt/UTIL_PION/REPORT_OUTPUT" ]; then
	echo "REPORT_OUTPUT sym link exits but is broken, replacing"
	rm "${USERPATH}/hallc_replay_lt/UTIL_PION/REPORT_OUTPUT"
	ln -s "${VOLATILEPATH}/${USER}/REPORT_OUTPUT/" "${USERPATH}/hallc_replay_lt/UTIL_PION/REPORT_OUTPUT"
    else echo "${USERPATH}/hallc_replay_lt/UTIL_PION/REPORT_OUTPUT sym link already exists and not broken"
    fi
fi

if [ ! -L "${USERPATH}/hallc_replay_lt/UTIL_PION/HISTOGRAMS" ]; then
    ln -s "${VOLATILEPATH}/${USER}/HISTOGRAMS/" "${USERPATH}/hallc_replay_lt/UTIL_PION/HISTOGRAMS"
elif [ -L "${USERPATH}/hallc_replay_lt/UTIL_PION/HISTOGRAMS" ]; then
    if [ ! -e "${USERPATH}/hallc_replay_lt/UTIL_PION/HISTOGRAMS" ]; then
	echo "HISTOGRAMS sym link exits but is broken, replacing"
	rm "${USERPATH}/hallc_replay_lt/UTIL_PION/HISTOGRAMS"
	ln -s "${VOLATILEPATH}/${USER}/HISTOGRAMS/" "${USERPATH}/hallc_replay_lt/UTIL_PION/HISTOGRAMS"
    else echo "${USERPATH}/hallc_replay_lt/UTIL_PION/HISTOGRAMS sym link already exists and not broken"
    fi
fi

echo ""
echo "Directories and sym links for UTIL_PION now setup"

exit 0
