#! /bin/bash
# Executes the python analysis script
#################################################################################################################################################

# source ROOT version 6.18.04
source /apps/root/6.18.04/setroot_CUE.bash

echo "Starting analysis of Pion events"
echo "I take as arguments the run prefix, run number and max number of events!"

RUNLIST=$1
if [[ -z "$1" ]]; then
    echo "I need a input RunList"
    echo "Please provide a run list as input"
fi

RUNPREFIX=$2
if [[ -z "$2" ]]; then
    echo "Please provide a run prefix as input"
    echo "Otherwise I'll assume coin_replay_Full"
    RUNPREFIX=coin_replay_Full
fi
MAXEVENTS=$3
if [[ -z "$3" ]]; then
    echo "Only Run Number entered...I'll assume -1 (all) events!" 
    MAXEVENTS=-1 
fi

#################################################################################################################################################

# Set path depending upon hostname. Change or add more as needed  
if [[ "${HOSTNAME}" = *"farm"* ]]; then  
    REPLAYPATH="/group/c-kaonlt/USERS/${USER}/hallc_replay_lt"
    cd "$REPLAYPATH"
    source "$REPLAYPATH/setup.sh"
elif [[ "${HOSTNAME}" = *"qcd"* ]]; then
    REPLAYPATH="/group/c-kaonlt/USERS/${USER}/hallc_replay_lt"
    cd "$REPLAYPATH"
    source "$REPLAYPATH/setup.sh" 
elif [[ "${HOSTNAME}" = *"cdaq"* ]]; then
    REPLAYPATH="/home/cdaq/hallc-online/hallc_replay_lt"
elif [[ "${HOSTNAME}" = *"phys.uregina.ca"* ]]; then
    REPLAYPATH="/home/${USER}/work/JLab/hallc_replay_lt"
fi
UTILPATH="${REPLAYPATH}/UTIL_PION"

###################################################################################################################################################

##Input run numbers##

inputFile="${UTILPATH}/scripts/work/${RUNLIST}"
while IFS='' read -r line || [[ -n "$line" ]]; do
    RUNNUMBER=$line
    if [ ! -f "${UTILPATH}/scripts/demo/OUTPUT/Analysis/PionLT/${RUNNUMBER}_${MAXEVENTS}_Analysed_Data.root" ]; then
        python3 ${UTILPATH}/scripts/work/analyzer.py ${RUNPREFIX} ${RUNNUMBER} ${MAXEVENTS}
    fi

done < "$inputFile"
