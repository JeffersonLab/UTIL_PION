#!/bin/bash

# 12/08/21 - Stephen JD Kay - University of Regina
# Script to fill the run list, will be executed by "master" script after each run

# Set up paths depending upon location

# Set path depending upon hostname. Change or add more as needed  
# Note, farm paths are only for teting
if [[ "${HOSTNAME}" = *"farm"* ]]; then  
    REPLAYPATH="/group/c-pionlt/USERS/${USER}/hallc_replay_lt"
elif [[ "${HOSTNAME}" = *"qcd"* ]]; then
    REPLAYPATH="/group/c-pionlt/USERS/${USER}/hallc_replay_lt"
elif [[ "${HOSTNAME}" = *"cdaq"* ]]; then
    REPLAYPATH="/home/cdaq/hallc-online/hallc_replay_lt"
fi

# Run number and run type should be read in by the "master" script, automating the target would be more difficult, master script should prompt for this - Only accept LD2, LH2, Dummy and carbon 1.5 or 6
RUNNUMBER=$1
RUNTYPE=$2
TARGET=$3
RUNLIST="${REPLAYPATH}/UTIL_PION/runlist_pionLT_2021.csv"
# Need to fix paths rather than give relative paths, also need to check information is still in these files and that it can grab it correctly
KINFILE="${REPLAYPATH}/DBASE/COIN/standard.kinematics"
# Need the names/paths for the files that contain the rest of the info we actually need
# SCALERFILE="OUTPUT/scalers_Run$RUNNUMBER.txt"
# REPORTFILE="../REPORT_OUTPUT/COIN/PRODUCTION/PionLT_replay_coin_production_${RUNNUMBER}_-1.report"
# MONITORFILE="../MON_OUTPUT/REPORT/reportMonitor_shms_${RUNNUMBER}_50000.txt"

# Get information available in standard.kinematics, execute a python script to do this for us
KINFILE_INFO=`python3 $REPLAYPATH/UTIL_PION/scripts/runlist/kinfile.py ${KINFILE} ${RUNNUMBER}`
# Split the string we get to individual variables, easier for printing and use later
SHMS_Angle=`echo ${KINFILE_INFO} | cut -d ','  -f1`
SHMS_P=`echo ${KINFILE_INFO} | cut -d ','  -f2`
HMS_Angle=`echo ${KINFILE_INFO} | cut -d ','  -f3`
HMS_P=`echo ${KINFILE_INFO} | cut -d ','  -f4`
EBeam=`echo ${KINFILE_INFO} | cut -d ','  -f5`

# Variables that still need to be set correctly
Current="Temp"
PS1="Temp"
PS4="Temp"
PS5="Temp"
HMS_Rate="Temp"
SHMS_Rate="Temp"
COIN_Rate="Temp"
Charge="Temp"
Raw_COIN="Temp"
Tracking="Temp"

 echo "========================================================================="
 echo "These values autofill into the run list ..."
 echo
 echo "Run number: $RUNNUMBER"
 echo "Run type: $RUNTYPE"
 echo "Target: $TARGET"
 echo "Beam energy: $EBeam"
 echo "SHMS momentum: $SHMS_P"
 echo "SHMS angle : $SHMS_Angle"
 echo "HMS momentum: $HMS_P"
 echo "HMS angle: $HMS_Angle"
# echo "Current: ${tmp[0]}"
# echo "PS1 : ${tmp[1]}"
# echo "PS3 : ${tmp[3]}"
# echo "PS5 : ${tmp[5]}"
# echo "HMS rate [kHz]: ${tmp[7]}"
# echo "SHMS rate [kHz]: ${tmp[8]}"
# echo "COIN rate [kHz]: ${tmp[9]}"
# echo "Charge [mC]: ${tmp[10]}"
# echo "Raw coin: ${tmp[11]}"
# echo "SHMS hadron tracking: ${tmp[12]}"
 echo "========================================================================="
# while true; do
#     read -p "Do these values all look correct? (Please answer yes or no) " yn
#     case $yn in
#         [Yy]* ) break;;
#         [Nn]* ) exit;;
#         * ) echo "Please answer yes or no.";;
#     esac
# done

# Ask user for a comment
read -p "Enter number of pi/n events and any other comments: " Comment
# Need to fix widths of entries with blank space at some point, see the test file for widths (based on headers)
RUNLIST_INFO="${RUNNUMBER},${RUNTYPE},${TARGET},${EBeam},${SHMS_P},${SHMS_Angle},${HMS_P},${HMS_Angle},${Current},${PS1},${PS4},${PS5},${HMS_Rate},${SHMS_Rate},${COIN_Rate},${Charge},${Raw_COIN},${Tracking},${Comment}"

# Check if there is already an entry for this run number, if there is, ask if you want to overwrite it, if not, print it to the file
DuplicateLines=() # Array to store line numbers of duplicated entries
LineNum=1 # Counter, starts at 1 since we skip the header
# Read run list, find any lines which already include an entry for this run number
while IFS='' read -r line || [[ -n "$line" ]]; do
    LineNum=$(($LineNum + 1 ))
    if [[ `echo ${line} | cut -d ','  -f1` == ${RUNNUMBER} ]]; then
	DuplicateLines[${#DuplicateLines[@]}]="${LineNum}"
    fi
done <  <(tail -n +2 ${RUNLIST}) # Ignores header line by using tail here

if [[ `echo "${#DuplicateLines[@]}"` != 0 ]]; then
    # Ask if the user wants to remove duplicate lines, in a grammatically correct manner :)
    if [[ `echo "${#DuplicateLines[@]}"` == 1 ]]; then
	read -p "$(echo "${#DuplicateLines[@]}") entry already found in the runlist for run ${RUNNUMBER}, delete dupliacte entry and print new entry to file? <Y/N> " prompt
    elif [[ `echo "${#DuplicateLines[@]}"` -gt 1 ]]; then
	read -p "$(echo "${#DuplicateLines[@]}") entries already found in the runlist for run ${RUNNUMBER}, delete dupliacte entries and print new entry to file? <Y/N> " prompt
    fi
    if [[ $prompt == "y" || $prompt == "Y" || $prompt == "yes" || $prompt == "Yes" ]]; then
	DeletedLines=0 # Counter to check how many lines we have deleted already
	# Loop over all line numbers identified earlier as being duplicates, delete them with a sed command
	for value in "${DuplicateLines[@]}"
	do
	    LineNum=$(($value-$DeletedLines)) # We need to account for any lines we delete as we go
	    sed -i "${LineNum}d" ${RUNLIST}
	    DeletedLines=$(($DeletedLines + 1))
	done
    echo ${RUNLIST_INFO} >> ${RUNLIST} # Print the run list info to the file
    else echo "Will not remove duplicate entries or print new entry to the file, please edit the runlist manually"
    fi
else
    echo ${RUNLIST_INFO} >> ${RUNLIST} # Print the run list info to the file
fi
