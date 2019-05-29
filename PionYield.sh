#!/bin/bash

echo "Starting Kaon Yield Estimation"
echo "I take as arguments the Run Number and max number of events!"
RUNNUMBER=$1
MAXEVENTS=$2
if [[ $1 -eq "" ]]; then
    echo "I need a Run Number!"
    exit 2
fi
if [[ $2 -eq "" ]]; then
    echo "Only Run Number entered...I'll assume -1 events!" 
    MAXEVENTS=-1 
fi
REPLAYPATH="/u/group/c-kaonlt/USERS/${USER}/hallc_replay_lt"
cd "$REPLAYPATH"
echo -e "\n\nStarting Replay Script\n\n"
eval "$REPLAYPATH/hcana -l -q \"UTIL_PION/scripts_Replay/replay_production_coin.C($RUNNUMBER,$MAXEVENTS)\"" | tee REPORT_OUTPUT/COIN/PRODUCTION/output_coin_production_${RUNNUMBER}_${MAXEVENTS}.report
cd "$REPLAYPATH/UTIL_PION/scripts_Yield/"
echo -e "\n\nYield Calculation\n\n"
root -l "run_PionYield.C($RUNNUMBER,$MAXEVENTS,5,1)"
