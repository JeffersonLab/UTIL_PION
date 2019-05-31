#!/bin/bash

echo "This script is for the analysis of Heep COINCIDENCE data!"
echo "Starting p(e,e'p) Yield Estimation"
echo "I take as arguments the Run Number and max number of events!"
RUNNUMBER=$1
MAXEVENTS=$2
if [[ $1 -eq "" ]]; then
    echo "I need a Run Number!"
    echo "Please provide a run number as input"
    exit 2
fi
if [[ $2 -eq "" ]]; then
    echo "Only Run Number entered...I'll assume -1 events!" 
    MAXEVENTS=-1 
fi

# Set path depending upon hostname. Change or add more as needed  
if [[ "${HOSTNAME}" = *"farm"* ]]; then  
    REPLAYPATH="/u/group/c-kaonlt/USERS/${USER}/hallc_replay_lt"
    if [[ "${HOSTNAME}" != *"ifarm"* ]]; then
	source /site/12gev_phys/softenv.sh 2.1
    fi
elif [[ "${HOSTNAME}" = *"cdaq"* ]]; then
    REPLAYPATH="/home/cdaq/hallc-online/hallc_replay_lt"
elif [[ "${HOSTNAME}" = *"phys.uregina.ca"* ]]; then
    REPLAYPATH="/home/${USER}/work/JLab/hallc_replay_lt"
fi
cd "/u/group/c-kaonlt/hcana/"
source "/u/group/c-kaonlt/hcana/setup.sh"
cd "$REPLAYPATH"
source "$REPLAYPATH/setup.sh"

eval "$REPLAYPATH/hcana -l -q \"UTIL_PION/scripts_Replay/replay_production_coin.C($RUNNUMBER,$MAXEVENTS)\""
sleep 15
cd "$REPLAYPATH/UTIL_PION/scripts_HeepCoin/"
if [[ "${HOSTNAME}" = *"farm"* && "${HOSTNAME}" != *"ifarm"* ]]; then
    root -l -b -q "run_HeepCoinYield.C($RUNNUMBER,$MAXEVENTS,5,1)"
else
    root -l "run_HeepCoinYield.C($RUNNUMBER,$MAXEVENTS,5,1)"
fi
exit 1
