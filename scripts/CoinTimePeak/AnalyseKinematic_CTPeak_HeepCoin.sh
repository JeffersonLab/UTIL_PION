#!/bin/bash
#19/10/20 - Stephen Kay, University of Regina

KINEMATIC=$1

if [[ -z "$1" ]]; then
    echo "I need a kinematic setting to process!"
    echo "Please provide a kinematic setting as input"
    exit 2
fi
declare -i Autosub=0
read -p "Auto submit batch jobs for missing replays/analyses? <Y/N> " prompt
if [[ $prompt == "y" || $prompt == "Y" || $prompt == "yes" || $prompt == "Yes" ]]; then
    Autosub=$((Autosub+1))
else echo "Will not submit any batch jobs, please check input lists manually and submit if needed"
fi

echo "######################################################"
echo "### Processing kinematic ${KINEMATIC} ###"
echo "######################################################"

# Set path depending upon hostname. Change or add more as needed  
if [[ "${HOSTNAME}" = *"farm"* ]]; then  
    REPLAYPATH="/group/c-pionlt/USERS/${USER}/hallc_replay_lt"
    if [[ "${HOSTNAME}" != *"ifarm"* ]]; then
	source /site/12gev_phys/softenv.sh 2.3
	source /apps/root/6.18.04/setroot_CUE.bash
    fi
    cd "$REPLAYPATH"
    source "$REPLAYPATH/setup.sh"
elif [[ "${HOSTNAME}" = *"qcd"* ]]; then
    REPLAYPATH="/group/c-pionlt/USERS/${USER}/hallc_replay_lt"
    source /site/12gev_phys/softenv.sh 2.3
    source /apps/root/6.18.04/setroot_CUE.bash
    cd "$REPLAYPATH"
    source "$REPLAYPATH/setup.sh" 
elif [[ "${HOSTNAME}" = *"cdaq"* ]]; then
    REPLAYPATH="/home/cdaq/hallc-online/hallc_replay_lt"
elif [[ "${HOSTNAME}" = *"phys.uregina.ca"* ]]; then
    REPLAYPATH="/home/${USER}/work/JLab/hallc_replay_lt"
fi
UTILPATH="${REPLAYPATH}/UTIL_PION"
RunListFile="${UTILPATH}/scripts/CoinTimePeak/Kinematics/${KINEMATIC}"
if [ ! -f "${RunListFile}" ]; then
    echo "Error, ${RunListFile} not found, exiting"
    exit 3
fi
cd $REPLAYPATH

if [ -f "${UTILPATH}/scripts/CoinTimePeak/Kinematics/${KINEMATIC}_MissingCTAnalysis" ]; then
    rm "${UTILPATH}/scripts/CoinTimePeak/Kinematics/${KINEMATIC}_MissingCTAnalysis"
else touch "${UTILPATH}/scripts/CoinTimePeak/Kinematics/${KINEMATIC}_MissingCTAnalysis"
fi

TestingVar=$((1))
while IFS='' read -r line || [[ -n "$line" ]]; do
    runNum=$line
    if [ ! -f "${UTILPATH}/OUTPUT/Analysis/PionLT/${runNum}_-1_CTPeak_Data_HeepCoin.root" ]; then
	echo "CTPeak analysis not found for run $runNum in ${UTILPATH}/scripts/CoinTimePeak/OUTPUT/"
	echo "${runNum}" >> "${UTILPATH}/scripts/CoinTimePeak/Kinematics/${KINEMATIC}_MissingCTAnalysis"
	TestingVar=$((TestingVar+1))
    fi
done < "$RunListFile"

if [ $TestingVar == 1 ]; then
    echo "All CT analysis files found"
    rm "${UTILPATH}/scripts/CoinTimePeak/Kinematics/${KINEMATIC}_MissingCTAnalysis"
elif [ $TestingVar != 1 ]; then
    cp "${UTILPATH}/scripts/CoinTimePeak/Kinematics/${KINEMATIC}_MissingCTAnalysis" "$REPLAYPATH/UTIL_BATCH/InputRunLists/${KINEMATIC}_MissingCTAnalysis"
    if [ $Autosub == 1 ]; then
	while IFS='' read -r line || [[ -n "$line" ]]; do
	    runNum=$line
	    if [ -f "${UTILPATH}/OUTPUT/Analysis/PionLT/${runNum}_-1_CTPeak_Data_HeepCoin.root" ]; then
		rm "${UTILPATH}/OUTPUT/Analysis/PionLT/${runNum}_-1_CTPeak_Data_HeepCoin.root"
	    fi
	    if [ -f "${UTILPATH}/ROOTfiles/Analysis/PionLT/Pion_coin_replay_production_${runNum}_-1.root" ]; then
		rm "${UTILPATH}/ROOTfiles/Analysis/PionLT/Pion_coin_replay_production_${runNum}_-1.root"
	    fi
	done < "${UTILPATH}/scripts/CoinTimePeak/Kinematics/${KINEMATIC}_MissingCTAnalysis"
	yes y | eval "$REPLAYPATH/UTIL_BATCH/batch_scripts/run_batch_CTPeak_HeepCoin_Analysis.sh ${KINEMATIC}_MissingCTAnalysis"
    elif [ $Autosub != 1 ]; then
	echo "Analyses missing, list copied to UTIL_BATCH directory, run on farm if desired"
	read -p "Process python script for missing replays/analyses interactively? <Y/N> " prompt2
	if [[ $prompt2 == "y" || $prompt2 == "Y" || $prompt2 == "yes" || $prompt2 == "Yes" ]]; then
	    while IFS='' read -r line || [[ -n "$line" ]]; do
		runNum=$line
		if [ ! -f "${UTILPATH}/OUTPUT/Analysis/PionLT/${runNum}_-1_CTPeak_Data_HeepCoin.root" ]; then
		    python3 $UTILPATH/scripts/CoinTimePeak/src/CoinTimePeak_HeepCoin.py "Pion_coin_replay_production" ${runNum} "-1" 
		fi
	    done < "$RunListFile"
	    else echo "Not processing python script interactively"
	fi
    fi
fi

if [ $TestingVar == 1 ]; then
    if [ -f "${UTILPATH}/OUTPUT/Analysis/PionLT/${KINEMATIC}_Output.csv" ]; then
	rm "${UTILPATH}/OUTPUT/Analysis/PionLT/${KINEMATIC}_Output.csv"
    else touch "${UTILPATH}/OUTPUT/Analysis/PionLT/${KINEMATIC}_Output.csv"
    fi
    while IFS='' read -r line || [[ -n "$line" ]]; do
	runNum=$line
	OutputFile="${UTILPATH}/OUTPUT/Analysis/PionLT/${runNum}_Out_tmp"
	if [ -f ${OutputFile} ]; then
	    rm ${OutputFile}
	else touch ${OutputFile}
	fi
	root -b -l -q "${UTILPATH}/scripts/CoinTimePeak/PlotHeepCoinPeak.C(\"${runNum}_-1_CTPeak_Data_HeepCoin.root\", \"${runNum}_CTOut_HeepCoin\")" >> ${OutputFile}
	sleep 1
	Data=$(sed -n "/${runNum},/p" $OutputFile)
	echo ${Data} >> "${UTILPATH}/OUTPUT/Analysis/PionLT/${KINEMATIC}_Output.csv"
	rm ${OutputFile}
	sleep 1
    done < "$RunListFile"
fi

exit 0

