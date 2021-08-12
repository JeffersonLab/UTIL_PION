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
# KINFILE="../DBASE/COIN/standard.kinematics"
# SCALERFILE="OUTPUT/scalers_Run$RUNNUMBER.txt"
# REPORTFILE="../REPORT_OUTPUT/COIN/PRODUCTION/PionLT_replay_coin_production_${RUNNUMBER}_-1.report"
# MONITORFILE="../MON_OUTPUT/REPORT/reportMonitor_shms_${RUNNUMBER}_50000.txt"
# SHMSANGLE=$(sed -n -e 's/^.*ptheta_lab = //p' $KINFILE | tail -1)
# SHMSMOMENT=$(sed -n -e 's/^.*ppcentral = //p' $KINFILE | tail -1)
# HMSANGLE=$(sed -n -e 's/^.*htheta_lab = //p' $KINFILE | tail -1)
# HMSMOMENT=$(sed -n -e 's/^.*hpcentral = //p' $KINFILE | tail -1)
# EBEAM=$(sed -n -e 's/^.*gpbeam = //p' $KINFILE | tail -1)
# python runlistGatherer.py $RUNNUMBER $RUNTYPE $TARGET
# inputFile="tmp"
# if [[ -e "output.txt" ]]; then
#     echo "Please close report file for replay and try again!"
#     exit 2
# fi
# if [[ $# -lt 3 ]]; then
#     echo "Not enough arguments, need run number, run type, and target"
#     exit 2
# elif [[ $# -gt 3 ]]; then
#     echo "Too many arguments, need only run number, run type, and target"
#     exit 2
# fi
# tmp=()
# while IFS='' read -r line || [[ -n "$line" ]]; do
#     tmp+=("$line")
# done < "$inputFile"
# echo "========================================================================="
# echo "These values autofill into the run list (i.e. emacs window)..."
# echo
# echo "Run number: $RUNNUMBER"
# echo "Run type: $RUNTYPE"
# echo "SHMS momentum: $SHMSMOMENT"
# echo "SHMS angle : $SHMSANGLE"
# echo "HMS momentum: $HMSMOMENT"
# echo "HMS angle: $HMSANGLE"
# echo "Target: $TARGET"
# echo "Beam energy: $EBEAM"
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
# echo "-------------------------------------------------------------------------"
# echo "You put these into the run sheet by hand (eek!)..."
# echo
# echo "SHMS prescale pretrig: ${tmp[13]}"
# echo "HMS prescale pretrig: ${tmp[14]}"
# echo "Coin pretrig: ${tmp[15]}"
# echo "========================================================================="
# while true; do
#     read -p "Do these values all look correct? (Please answer yes or no) " yn
#     case $yn in
#         [Yy]* ) break;;
#         [Nn]* ) exit;;
#         * ) echo "Please answer yes or no.";;
#     esac
# done
# read -p "Enter neutron events and any other comments: " comment
# echo -e "$RUNNUMBER\t$RUNTYPE\t$SHMSMOMENT\t$SHMSANGLE\t\t$HMSMOMENT\t$HMSANGLE\t\t$TARGET\t$EBEAM\t${tmp[0]}\t${tmp[1]}\t${tmp[3]}\t${tmp[5]}\t${tmp[7]}\t\t${tmp[8]}\t\t${tmp[9]}\t\t${tmp[10]}\t\t${tmp[11]}\t\t${tmp[12]}\t\t\t!$comment" >> $RUNLIST
    
