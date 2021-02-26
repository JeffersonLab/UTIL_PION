#!/bin/bash

# 07/05/20 - Stephen Kay, University of Regina
# Runs the analysis for the setting indicated, in this case it is split into two halves due to a change in the beam conditions
RUNPREFIX=$1
if [[ -z "$1" ]]; then
    echo "I need a Run Prefix!"
    echo "Please provide a run prefix as input"
    exit 2
fi

# Simply hadd the two halves together if you want

echo "Starting analysis of Q2 = 4.4, W = 2.74, central angle, low espilon setting"

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
RunListFile="${UTILPATH}/scripts/kinematics/Q4p4W2p74center_lowe_Pt1"
while IFS='' read -r line || [[ -n "$line" ]]; do
    runNum=$line
    RootName+="${runNum}_-1_Analysed_Data.root "
    eval '"$SCRIPTPATH" $RunPrefix $runNum -1'
done < "$RunListFile"
sleep 5
cd "${UTILPATH}/OUTPUT/Analysis/PionLT"
KINFILE="Q4p4W2p74center_lowe_Pt1.root"
hadd ${KINFILE} ${RootName}

if [ ! -f "${UTILPATH}/OUTPUT/Analysis/PionLT/Q4p4W2p74center_lowe_Pt1_Pions.root" ]; then
    root -b -l -q "${UTILPATH}/scripts/pionyield/PlotPionPhysics.C(\"${KINFILE}\", \"Q4p4W2p74center_lowe_Pt1_Kaons\")"
elif [ ! -f "${UTILPATH}/OUTPUT/Analysis/PionLT/Q4p4W2p74center_lowe_Pt1_Pions.pdf" ]; then
    root -b -l -q "${UTILPATH}/scripts/pionyield/PlotPionPhysics.C(\"${KINFILE}\", \"Q4p4W2p74center_lowe_Pt1_Kaons\")"
else echo "Pion plots already found in - ${UTILPATH}/OUTPUT/Analysis/PionLT/Q4p4W2p74center_lowe_Pt1_Kaons.root and .pdf - Plotting macro skipped"
fi
exit 0
