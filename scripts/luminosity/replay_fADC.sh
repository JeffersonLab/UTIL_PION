#!/bin/bash

echo "Starting fADC Script"
echo "I take as arguments the Run Number and max number of events!"
RUNNUMBER=$1
MAXEVENTS=$2
#MAXEVENTS=12500

if [[ $1 -eq "" ]]; then
    echo "I need a Run Number!"
    exit 2
fi

if [[ $2 -eq "" ]]; then
    echo "Only Run Number entered...I'll assume -1 events!" 
    MAXEVENTS=-1 
fi

# Set path depending upon hostname. Change or add more as needed  
if [[ "${HOSTNAME}" = *"farm"* ]]; then  
    #REPLAYPATH="/group/c-pionlt/USERS/${USER}/hallc_replay_lt"    
    REPLAYPATH="/group/c-pionlt/online_analysis/hallc_replay_lt"
    if [[ "${HOSTNAME}" != *"ifarm"* ]]; then
	source /site/12gev_phys/softenv.sh 2.3
	source /apps/root/6.18.04/setroot_CUE.bash
    fi
    cd "/group/c-pionlt/hcana/"
    source "/group/c-pionlt/hcana/setup.sh"
    cd "$REPLAYPATH"
    source "$REPLAYPATH/setup.sh"
elif [[ "${HOSTNAME}" = *"qcd"* ]]; then
    #REPLAYPATH="/group/c-pionlt/USERS/${USER}/hallc_replay_lt"
    REPLAYPATH="/group/c-pionlt/online_analysis/hallc_replay_lt"
    source /site/12gev_phys/softenv.sh 2.3
    source /apps/root/6.18.04/setroot_CUE.bash
    cd "/group/c-pionlt/hcana/"
    source "/group/c-pionlt/hcana/setup.sh" 
    cd "$REPLAYPATH"
    source "$REPLAYPATH/setup.sh" 
elif [[ "${HOSTNAME}" = *"cdaq"* ]]; then
    REPLAYPATH="/home/cdaq/hallc-online/hallc_replay_lt"
elif [[ "${HOSTNAME}" = *"phys.uregina.ca"* ]]; then
    REPLAYPATH="/home/${USER}/work/JLab/hallc_replay_lt"
elif [[ "${HOSTNAME}" = *"trottar"* ]]; then
    REPLAYPATH="/home/trottar/Analysis/hallc_replay_lt"
fi

UTILPATH="${REPLAYPATH}/UTIL_PION"
cd $REPLAYPATH

###################################################################################################################################################
# RLT 09/24/21...Changed from 150k to full analysis. There may be issues with the currents/EDTM because the cuts may only be applying trip cuts to the first 150k events.
# Section for luminosity replay script
if [ ! -f "$REPLAYPATH/UTIL_PION/ROOTfiles/Scalers/coin_replay_scalers_${RUNNUMBER}_${MAXEVENTS}.root" ]; then
    eval "$REPLAYPATH/hcana -l -q -b \"SCRIPTS/COIN/SCALERS/replay_coin_scalers.C($RUNNUMBER,${MAXEVENTS})\""
    cd "$REPLAYPATH/CALIBRATION/bcm_current_map"
    root -b -l<<EOF 
.L ScalerCalib.C+
.x run.C("${REPLAYPATH}/UTIL_PION/ROOTfiles/Scalers/coin_replay_scalers_${RUNNUMBER}_${MAXEVENTS}.root")
.q  
EOF
    mv bcmcurrent_$RUNNUMBER.param $REPLAYPATH/PARAM/HMS/BCM/CALIB/bcmcurrent_$RUNNUMBER.param
    cd $REPLAYPATH
else echo "Scaler replayfile already found for this run in $REPLAYPATH/ROOTfiles/Scalers - Skipping scaler replay step"
fi

sleep 3
# SJDK 31/08/21 - Replays for luminosity analysis should output to Analysis/Lumi, for now this is probably fine
if [ ! -f "$REPLAYPATH/UTIL_PION/ROOTfiles/Analysis/Lumi/Pion_replay_luminosity_${RUNNUMBER}_${MAXEVENTS}.root" ]; then
    if [[ "${HOSTNAME}" != *"ifarm"* ]]; then
	if [[ "${HOSTNAME}" == *"cdaq"* ]]; then
	    eval "$REPLAYPATH/hcana -l -q -b \"UTIL_PION/scripts/replay/replay_luminosity.C($RUNNUMBER,$MAXEVENTS)\""| tee $REPLAYPATH/UTIL_PION/REPORT_OUTPUT/Analysis/Lumi/Pion_output_coin_production_Summary_${RUNNUMBER}_${MAXEVENTS}.report
	else	
	    eval "$REPLAYPATH/hcana -l -q -b \"UTIL_PION/scripts/replay/replay_luminosity.C($RUNNUMBER,$MAXEVENTS)\"" 
	fi
    elif [[ "${HOSTNAME}" == *"ifarm"* ]]; then
	eval "$REPLAYPATH/hcana -l -q -b \"UTIL_PION/scripts/replay/replay_luminosity.C($RUNNUMBER,$MAXEVENTS)\""| tee $REPLAYPATH/UTIL_PION/REPORT_OUTPUT/Analysis/Lumi/Pion_output_coin_production_Summary_${RUNNUMBER}_${MAXEVENTS}.report
    fi
else echo "Replayfile already found for this run in $REPLAYPATH/UTIL_PION/ROOTfiles/Analysis/Lumi/ - Skipping replay step"
fi

sleep 3

source /site/12gev_phys/softenv.sh 2.3
source /apps/root/6.18.04/setroot_CUE.bash

# Sets trigger windows
echo
echo "Running trigWindows.sh ${RUNNUMBER}..."
cd ${REPLAYPATH}/UTIL_PION/scripts/trig_windows/
source trigWindows.sh ${RUNNUMBER}

# Analyzes lumi runs
echo "Running lumiyield.py ${RUNNUMBER} ${MAXEVENTS}..."
cd ${REPLAYPATH}/UTIL_PION/scripts/luminosity/src/
python3 fADCyield.py Pion_replay_luminosity ${RUNNUMBER} ${MAXEVENTS}
