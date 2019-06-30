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
    REPLAYPATH="/group/c-kaonlt/USERS/${USER}/hallc_replay_lt"
    if [[ "${HOSTNAME}" != *"ifarm"* ]]; then
	source /site/12gev_phys/softenv.sh 2.1
    fi
    cd "/group/c-kaonlt/hcana/"
    source "/group/c-kaonlt/hcana/setup.sh"
    cd "$REPLAYPATH"
    source "$REPLAYPATH/setup.sh"
elif [[ "${HOSTNAME}" = *"cdaq"* ]]; then
    REPLAYPATH="/home/cdaq/hallc-online/hallc_replay_lt"
elif [[ "${HOSTNAME}" = *"phys.uregina.ca"* ]]; then
    REPLAYPATH="/home/${USER}/work/JLab/hallc_replay_lt"
fi
cd $REPLAYPATH
if [ ! -f "$REPLAYPATH/ROOTfilesPion/PionLT_coin_replay_production_${RUNNUMBER}_${MAXEVENTS}.root" ]; then
    eval "$REPLAYPATH/hcana -l -q \"UTIL_PION/scripts_Replay/replay_production_coin.C($RUNNUMBER,$MAXEVENTS)\"" | tee $REPLAYPATH/UTIL_PION/REPORT_OUTPUT/COIN/PRODUCTION/PionLT_output_coin_production_${RUNNUMBER}_${MAXEVENTS}.report 
fi
sleep 5
cd "$REPLAYPATH/UTIL_PION/scripts_HeepCoin/"
if [[ "${HOSTNAME}" = *"farm"* && "${HOSTNAME}" != *"ifarm"* ]]; then
    root -l -b -q "run_HeepCoinYield.C($RUNNUMBER,$MAXEVENTS,5,1)"
else
    root -l "run_HeepCoinYield.C($RUNNUMBER,$MAXEVENTS,5,1)"
fi
if [[ "${HOSTNAME}" = *"cdaq"* ]]; then
    cd "$REPLAYPATH/UTIL_PION"
    python reportSummary.py $RUNNUMBER $MAXEVENTS
    emacs output.txt
    mv output.txt OUTPUT/scalers_Run$RUNNUMBER.txt
    if [[ -e "OUTPUT/scalers_Run$RUNNUMBER.txt" ]]; then
	while true; do
	    read -p "Would you like to update the run list as well? (Please answer yes or no) " yn
	    case $yn in
		[Yy]* ) break;;
		[Nn]* ) exit;;
		* ) echo "Please answer yes or no.";;
	    esac
	done
	read -p "What type of production was this run? (e.g. Prod, Heep, ect.)" runType
	read -p "What was the target? (Dummy, LH2 ..?)" target
	fillrunList="./fill_runList $RUNNUMBER $runType $target"
	eval ${fillrunList}
    fi
fi
exit 0
