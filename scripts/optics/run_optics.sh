#! /bin/bash

# Stephen JD Kay - University of Regina - 07/09/21
# This script executes a full (all events) version of the 50k replay on each spectrometer arm
# It then produces a few quick plots of the relevant optics quantities
# To run this script, execute ./run_optics $RUNNUMBER$

#################################################################################################################################################

echo "Starting optics analysis"
echo "I take 1 argument - run number"
# Input params - run number
RUNNUMBER=$1
if [[ -z "$1" ]]; then
    echo "I need an input run number"
    echo "Please provide a run number as input"
fi

#################################################################################################################################################

# Set path depending upon hostname. Change or add more as needed  
if [[ "${HOSTNAME}" = *"cdaq"* ]]; then
    REPLAYPATH="/home/cdaq/hallc-online/hallc_replay_lt"
fi

UTILPATH="${REPLAYPATH}/UTIL_PION"
cd $REPLAYPATH

###################################################################################################################################################

# Carry out a full replay of each spectrometer
if [ ! -f "$REPLAYPATH/ROOTfiles/Analysis/50k/hms_coin_replay_production_${RUNNUMBER}_-1.root" ]; then
    eval "$REPLAYPATH/hcana -l -q \"SCRIPTS/HMS/PRODUCTION/replay_production_hms_coin.C($RUNNUMBER,-1)\""
else echo "HMS replayfile already found for this run in $REPLAYPATH/ROOTfiles/Analysis/50k/ - Skipping HMS replay step"
fi

sleep 3

if [ ! -f "$REPLAYPATH/ROOTfiles/Analysis/50k/shms_coin_replay_production_${RUNNUMBER}_-1.root" ]; then
    eval "$REPLAYPATH/hcana -l -q \"SCRIPTS/SHMS/PRODUCTION/replay_production_shms_coin.C($RUNNUMBER,-1)\""
else echo "SHMS replayfile already found for this run in $REPLAYPATH/ROOTfiles/Analysis/50k/ - Skipping SHMS replay step"
fi

exit 0
