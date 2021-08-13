#! /bin/bash

# 10/July/2021, Author - Muhammad Junaid, University of Regina, Canada

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

#RUNPREFIX=$2
#if [[ -z "$2" ]]; then
#    echo "Please provide a run prefix as input"
#    echo "Otherwise I'll assume coin_replay_Full"
#fi

MAXEVENTS=$2
if [[ -z "$2" ]]; then
    echo "Only Run Number entered...I'll assume -1 (all) events!" 
    MAXEVENTS=-1 
fi

#################################################################################################################################################

# Set path depending upon hostname. Change or add more as needed  
if [[ "${HOSTNAME}" = *"farm"* ]]; then  
    REPLAYPATH="/group/c-kaonlt/USERS/${USER}/hallc_replay_lt"
    if [[ "${HOSTNAME}" != *"ifarm"* ]]; then
	source /site/12gev_phys/softenv.sh 2.1
    fi
    cd "/group/c-kaonlt/hcana/"
    source "/group/c-kaonlt/hcana/setup.sh"
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
    if [ ! -f "${REPLAYPATH}/ROOTfiles/Analysis/General/coin_replay_Full_${RUNNUMBER}_${MAXEVENTS}.root" ]; then
        eval "${REPLAYPATH}/hcana -l -q \"SCRIPTS/COIN/PRODUCTION/FullReplay.C ($RUNNUMBER,$MAXEVENTS)\""
    else echo "Replay root file already found in ${REPLAYPATH}/ROOTfiles/Analysis/General/ - Skipped python replay script step"
    fi

    RUNNUMBER=$line
    if [ ! -f "${UTILPATH}/OUTPUT/Analysis/PionLT/${RUNNUMBER}_${MAXEVENTS}_Analysed_Data.root" ]; then
        python3 ${UTILPATH}/scripts/work/analyzer.py coin_replay_Full ${RUNNUMBER} ${MAXEVENTS}
    else echo "Analysed root file already found in ${UTILPATH}/OUTPUT/Analysis/PionLT/ - Skipped python analyzer script step"
    fi

    RUNNUMBER=$line
    if [ ! -f "${UTILPATH}/OUTPUT/Analysis/PionLT/${RUNNUMBER}_${MAXEVENTS}_Output_Data.root" ]; then
        python3 ${UTILPATH}/scripts/work/PlotPionPhysics.py ${RUNNUMBER} ${MAXEVENTS} Analysed_Data
    else echo "Output root file already found in ${UTILPATH}/OUTPUT/Analysis/PionLT/ - Skipped python output script step"
    fi

done < "$inputFile"
