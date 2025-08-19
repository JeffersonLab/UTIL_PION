#! /bin/bash

#
# Description:
# ================================================================
# Time-stamp: "2024-03-11 12:12:00 junaid"
# ================================================================
#
# Author:  Muhammad Junaid <mjo147@uregina.ca>
#
# Copyright (c) junaid
#
##################################################################################
# Created - 10/July/2021, Author - Muhammad Junaid, University of Regina, Canada
##################################################################################

echo "Starting analysis of Proton events"
echo "I take as arguments the beam energy and number of events!"
# Input params - beam energy and max number of events
BEAM_ENERGY=$1
if [[ -z "$1" ]]; then
    echo "I need a beam energy"
    echo "Please provide a beam energy as input"
fi
MAXEVENTS=$2
if [[ -z "$2" ]]; then
    echo "Only beam energy entered...I'll assume -1 (all) events!" 
    MAXEVENTS=-1 
fi
DATA_Suffix=HeePCoin_Analysed_Data
DUMMY_Suffix=HeePCoin_Analysed_Dummy_Data
SIMC_Suffix=Heep_Coin_SIMC
DATA_RUN_LIST=HeePCoin_${BEAM_ENERGY}
DUMMY_RUN_LIST=HeePCoin_${BEAM_ENERGY}_Dummy
CSV_FILE=PionLT_HeeP_coin_HeePCoin_efficiency_data_2025_03_08

# Runs script in the ltsep python package that grabs current path enviroment
if [[ ${HOSTNAME} = *"cdaq"* ]]; then
    PATHFILE_INFO=`python3 /home/cdaq/pionLT-2021/hallc_replay_lt/UTIL_PION/bin/python/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
elif [[ "${HOSTNAME}" = *"farm"* ]]; then
#    PATHFILE_INFO=`python3 /u/home/${USER}/.local/lib/python3.4/site-packages/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
    PATHFILE_INFO=`python3 /u/group/c-pionlt/USERS/junaid/replay_lt_env/lib/python3.9/site-packages/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
elif [[ "${HOSTNAME}" = *"qcd"* ]]; then
    PATHFILE_INFO=`python3 /u/home/${USER}/.local/lib/python3.4/site-packages/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
fi

# Split the string we get to individual variables, easier for printing and use later
VOLATILEPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f1` # Cut the string on , delimitter, select field (f) 1, set variable to output of command
ANALYSISPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f2`
HCANAPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f3`
REPLAYPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f4`
UTILPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f5`
PACKAGEPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f6`
OUTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f7`
ROOTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f8`
REPORTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f9`
CUTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f10`
PARAMPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f11`
SCRIPTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f12`
ANATYPE=`echo ${PATHFILE_INFO} | cut -d ','  -f13`
USER=`echo ${PATHFILE_INFO} | cut -d ','  -f14`
HOST=`echo ${PATHFILE_INFO} | cut -d ','  -f15`

#################################################################################################################################################

cd $REPLAYPATH

# Section for HeeP physics ploting script
if [ -f "${UTILPATH}/OUTPUT/Analysis/HeeP/${BEAM_ENERGY}_${MAXEVENTS}_Output_Data.root" ]; then
    read -p "HeeP coin output plots file already exists, do you want to reprocess it? <Y/N> " option2
    if [[ $option2 == "y" || $option2 == "Y" || $option2 == "yes" || $option2 == "Yes" ]]; then
	rm "${UTILPATH}/OUTPUT/Analysis/HeeP/${BEAM_ENERGY}_${MAXEVENTS}_Output_Data.root"
	echo "Reprocessing"
	python3 ${UTILPATH}/scripts/heep/src/plot_coin_comp.py ${BEAM_ENERGY} ${MAXEVENTS} ${DATA_Suffix} ${DUMMY_Suffix} ${SIMC_Suffix} ${DATA_RUN_LIST} ${DUMMY_RUN_LIST} ${CSV_FILE}
    else
	echo "Skipping python HeeP plotting script step"
    fi
elif [ ! -f  "${UTILPATH}/OUTPUT/Analysis/HeeP/${BEAM_ENERGY}_${MAXEVENTS}_Output_Data.root" ]; then
	python3 ${UTILPATH}/scripts/heep/src/plot_coin_comp.py ${BEAM_ENERGY} ${MAXEVENTS} ${DATA_Suffix} ${DUMMY_Suffix} ${SIMC_Suffix} ${DATA_RUN_LIST} ${DUMMY_RUN_LIST} ${CSV_FILE}
else echo "HeeP coin output plots file already found in ${UTILPATH}/OUTPUT/Analysis/HeeP/ - Skipped python plotting script step"
fi
#evince "${UTILPATH}/OUTPUT/Analysis/HeeP/${BEAM_ENERGY}_${MAXEVENTS}_heep_Proton_Analysis_Distributions.pdf" &
exit 0
