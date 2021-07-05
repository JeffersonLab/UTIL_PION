#!/bin/bash

# 26/05/20 - Stephen Kay, University of Regina

RUNPREFIX=$1
if [[ -z "$1" ]]; then
    echo "I need a Run Prefix!"
    echo "Please provide a run prefix as input"
    exit 2
fi

echo "Starting analysis of Q2 = 3.0, W = 3.14, central angle, high espilon setting"

# Set path depending upon hostname. Change or add more as needed  
if [[ "${HOSTNAME}" = *"farm"* ]]; then  
    REPLAYPATH="/group/c-pionlt/USERS/${USER}/hallc_replay_lt"
    if [[ "${HOSTNAME}" != *"ifarm"* ]]; then
	source /site/12gev_phys/softenv.sh 2.3
    fi
    cd "$REPLAYPATH"
    source "$REPLAYPATH/setup.sh"
elif [[ "${HOSTNAME}" = *"qcd"* ]]; then
    REPLAYPATH="/group/c-pionlt/USERS/${USER}/hallc_replay_lt"
    source /site/12gev_phys/softenv.sh 2.3
    cd "$REPLAYPATH"
    source "$REPLAYPATH/setup.sh" 
elif [[ "${HOSTNAME}" = *"cdaq"* ]]; then
    REPLAYPATH="/home/cdaq/hallc-online/hallc_replay_lt"
elif [[ "${HOSTNAME}" = *"phys.uregina.ca"* ]]; then
    REPLAYPATH="/home/${USER}/work/JLab/hallc_replay_lt"
fi
UTILPATH="${REPLAYPATH}/UTIL_PION"
SCRIPTPATH="${REPLAYPATH}/UTIL_PION/scripts/pionyield/Analyse_Pions.sh"
RunListFile="${UTILPATH}/scripts/kinematics/Q3W3p14center_highe"
while IFS='' read -r line || [[ -n "$line" ]]; do
    runNum=$line
    RootName+="${runNum}_-1_Analysed_Data.root "
    eval '"$SCRIPTPATH" $RunPrefix $runNum -1'
done < "$RunListFile"
sleep 5
cd "${UTILPATH}/OUTPUT/Analysis/PionLT"
KINFILE="Q3W3p14center_highe.root"
hadd ${KINFILE} ${RootName}

if [ ! -f "${UTILPATH}/OUTPUT/Analysis/PionLT/Q3W3p14center_highe_Pions.root" ]; then
    root -b -l -q "${UTILPATH}/scripts/pionyield/PlotPionPhysics.C(\"${KINFILE}\", \"Q3W3p14center_highe_Kaons\")"
elif [ ! -f "${UTILPATH}/OUTPUT/Analysis/PionLT/Q3W3p14center_highe_Pions.pdf" ]; then
    root -b -l -q "${UTILPATH}/scripts/pionyield/PlotPionPhysics.C(\"${KINFILE}\", \"Q3W3p14center_highe_Kaons\")"
else echo "Pion plots already found in - ${UTILPATH}/OUTPUT/Analysis/PionLT/Q3W3p14center_highe_Kaons.root and .pdf - Plotting macro skipped"
fi
exit 0
