#!/bin/bash

# 26/05/20 - Stephen Kay, University of Regina

RUNPREFIX=$1
if [[ -z "$1" ]]; then
    echo "I need a Run Prefix!"
    echo "Please provide a run prefix as input"
    exit 2
fi

echo "Starting analysis of Q2 = 4.4, W = 2.74, left angle, low espilon setting"

# Set path depending upon hostname. Change or add more as needed  
if [[ "${HOSTNAME}" = *"farm"* ]]; then  
    REPLAYPATH="/group/c-kaonlt/USERS/${USER}/hallc_replay_lt"
    if [[ "${HOSTNAME}" != *"ifarm"* ]]; then
	source /site/12gev_phys/softenv.sh 2.3
    fi
    cd "$REPLAYPATH"
    source "$REPLAYPATH/setup.sh"
elif [[ "${HOSTNAME}" = *"qcd"* ]]; then
    REPLAYPATH="/group/c-kaonlt/USERS/${USER}/hallc_replay_lt"
    source /site/12gev_phys/softenv.sh 2.3
    cd "$REPLAYPATH"
    source "$REPLAYPATH/setup.sh" 
elif [[ "${HOSTNAME}" = *"cdaq"* ]]; then
    REPLAYPATH="/home/cdaq/hallc-online/hallc_replay_lt"
elif [[ "${HOSTNAME}" = *"phys.uregina.ca"* ]]; then
    REPLAYPATH="/home/${USER}/work/JLab/hallc_replay_lt"
fi
UTILPATH="${REPLAYPATH}/UTIL_KAONLT"
SCRIPTPATH="${REPLAYPATH}/UTIL_KAONLT/scripts/kaonyield/Analyse_Kaons.sh"
RunListFile="${UTILPATH}/scripts/kinematics/Q4p4W2p74left_lowe"
while IFS='' read -r line || [[ -n "$line" ]]; do
    runNum=$line
    RootName+="${runNum}_-1_Analysed_Data.root "
    eval '"$SCRIPTPATH" $RunPrefix $runNum -1'
done < "$RunListFile"
sleep 5
cd "${UTILPATH}/scripts/kaonyield/OUTPUT"
KINFILE="Q4p4W2p74left_lowe.root"
hadd ${KINFILE} ${RootName}

if [ ! -f "${UTILPATH}/scripts/kaonyield/OUTPUT/Q4p4W2p74left_lowe_Kaons.root" ]; then
    root -b -l -q "${UTILPATH}/scripts/kaonyield/PlotKaonPhysics.C(\"${KINFILE}\", \"Q4p4W2p74left_lowe_Kaons\")"
elif [ ! -f "${UTILPATH}/scripts/kaonyield/OUTPUT/Q4p4W2p74left_lowe_Kaons.pdf" ]; then
    root -b -l -q "${UTILPATH}/scripts/kaonyield/PlotKaonPhysics.C(\"${KINFILE}\", \"Q4p4W2p74left_lowe_Kaons\")"
else echo "Kaon plots already found in - ${UTILPATH}/scripts/kaonyield/OUTPUT/Q4p4W2p74left_lowe_Kaons.root and .pdf - Plotting macro skipped"
fi
exit 0
