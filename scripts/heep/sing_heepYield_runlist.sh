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

while getopts 'cpb' flag; do
    case "${flag}" in
        h)
        echo "-------------------------------------------------------------------"
        echo "./coin_heepYield_runlist.sh -{flags} {variable arguments, see help}"
        echo "-------------------------------------------------------------------"
        echo
        echo "The following flags can be called for the heep analysis..."         
        echo "    -h, help"
        echo "    -c, run cut analyser script to trim down the root file for further analysis"
        echo "        cut -> RunType=arg1 Spec=arg2 RunList=arg3 MaxEvents=arg4"        
        echo "    -p, run plotting script to save output as pdf"
        echo "        plot -> RunType=arg1 Spec=arg2 RunList=arg3 MaxEvents=arg4"  
        echo "    -b, run cut analyser and plotting script to save output as pdf for beamenergy"
        echo "        plot -> RunType=arg1 Spec=arg2 RunList=arg3 MaxEvents=arg4"  
        exit 0
        ;;
        c) c_flag='true' ;;
        p) p_flag='true' ;;
        b) b_flag='true' ;;
        *) print_usage
        exit 1 ;;
    esac
done

echo "Starting analysis of Proton events"
echo "I take as arguments the -flag, beam energy, runlist and number of events!"

# Input params - beam energy, run list and max number of events
RunType=$2
if [[ -z "$2" ]]; then
    echo "No RunType provided. Please provide type of Analysis 'HeePSing'!"
fi

Spec=$3
if [[ -z "$3" ]]; then
    echo "No Spectrometer provided. Please provide SPEC type 'HMS' or 'SHMS'!"
fi

RunList=$4
if [[ -z "$4" ]]; then
    echo "I need a runlist"
    echo "Please provide a run list as input"
fi

MAXEVENTS=$5
if [[ -z "$5" ]]; then
    echo "Only Run list entered...I'll assume -1 (all) events!" 
    MAXEVENTS=-1 
fi

# Runs script in the ltsep python package that grabs current path enviroment
if [[ ${HOSTNAME} = *"cdaq"* ]]; then
    PATHFILE_INFO=`python3 /home/cdaq/pionLT-2021/hallc_replay_lt/UTIL_PION/bin/python/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
elif [[ "${HOSTNAME}" = *"farm"* ]]; then
#    PATHFILE_INFO=`python3 /u/home/${USER}/.local/lib/python3.4/site-packages/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
    PATHFILE_INFO=`python3 /u/group/c-pionlt/USERS/${USER}/replay_lt_env/lib/python3.9/site-packages/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
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

cd $REPLAYPATH
inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/heep_runlist/${RunList}"
BEAM_ENERGY=$(echo "${RunList}" | awk -F'_' '{print $2}')
DATA_Suffix=HeePSing_${Spec}_Analysed_Data
DUMMY_Suffix=HeePSing_${Spec}_Analysed_Dummy_Data
SIMC_Suffix=Heep_Sing_${Spec}_SIMC
DATA_RUN_LIST=HeePSing_${BEAM_ENERGY}
DUMMY_RUN_LIST=HeePSing_${BEAM_ENERGY}_Dummy
CSV_FILE=PionLT_${Spec}_HeePSing_HeePSing_efficiency_data_2024_03_02

#if [[ "$Spec" == "HMS" ]]; then
#     CSV_FILE=PionLT_${Spec}_HeePSing_hms_HeePCoin_efficiency_data_2024_09_16
#elif [[ "$Spec" == "SHMS" ]]; then
#     CSV_FILE=PionLT_HeePSing_shms_HeePCoin_efficiency_data_2024_09_16
#fi
################################################################################################################################                                                                                   

# Section for HeeP analysis script
ROOTPREFIX_CUT=${ANATYPE}LT_${Spec}_HeePSing_replay_production

if [[ $c_flag = "true" ]]; then
    while true; do
        read -p "Do you wish to analyse root files with runlist ${inputFile}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                i=-1
                (
                ##Reads in input file##
                while IFS='' read -r line || [[ -n "$line" ]]; do
                    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                    echo "Run number read from file: $line"
                    echo ""
                    cd "${UTILPATH}/scripts/heep/src/"
                    python3 singyield_heep.py ${ROOTPREFIX_CUT} $line ${MAXEVENTS} ${Spec}
                done < "$inputFile"
                cd $REPLAYPATH/OUTPUT/Analysis/HeeP/
                dir_name="runbyrun"
                # Check if the directory exists
                if [ ! -d "$dir_name" ]; then
                   # If it doesn't exist, create it
                   mkdir "$dir_name"
                   echo "Directory '$dir_name' created."
                else
                   echo "Directory '$dir_name' already exists."
                fi
#                if [[ "$RunList" == *"Dummy"* ]]; then
                if [[ "${RunList,,}" == *"dummy"* ]]; then
                    if [ -e "${BEAM_ENERGY}_${MAXEVENTS}_HeePSing_${Spec}_Analysed_Dummy_Data.root" ]; then
                        echo "Deleting ${BEAM_ENERGY}_${MAXEVENTS}_HeePSing_${Spec}_Analysed_Dummy_Data.root"
    		        rm ${BEAM_ENERGY}_${MAXEVENTS}_HeePSing_${Spec}_Analysed_Dummy_Data.root
                    fi  
		    hadd ${BEAM_ENERGY}_${MAXEVENTS}_HeePSing_${Spec}_Analysed_Dummy_Data.root 1*
		else
                    if [ -e "${BEAM_ENERGY}_${MAXEVENTS}_HeePSing_${Spec}_Analysed_Data.root" ]; then
			echo "Deleting ${BEAM_ENERGY}_${MAXEVENTS}_HeePSing_${Spec}_Analysed_Data.root"
			rm ${BEAM_ENERGY}_${MAXEVENTS}_HeePSing_${Spec}_Analysed_Data.root
                    fi
		    hadd ${BEAM_ENERGY}_${MAXEVENTS}_HeePSing_${Spec}_Analysed_Data.root 1*
         	fi
		mv 1* $dir_name
          	)
                break;;
            [Nn]* )
                exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

sleep 3

#################################################################################################################################
# Section for HeeP physics ploting script

elif [[ $p_flag = "true" ]]; then
    while true; do
        read -p "Do you wish to plot the ROOT files with runlist ${inputFile}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                i=-1
                (
                # Section for HeeP physics ploting script
                if [ -f "${UTILPATH}/OUTPUT/Analysis/HeeP/${BEAM_ENERGY}_${MAXEVENTS}_HeePSing_${Spec}_Output_Data.root" ]; then
                    read -p "HeeP coin output plots file already exists, do you want to reprocess it? <Y/N> " option2
                    if [[ $option2 == "y" || $option2 == "Y" || $option2 == "yes" || $option2 == "Yes" ]]; then
                        rm "${UTILPATH}/OUTPUT/Analysis/HeeP/${BEAM_ENERGY}_${MAXEVENTS}_HeePSing_${Spec}_Output_Data.root"
                        echo "Reprocessing"
                        python3 ${UTILPATH}/scripts/heep/src/plot_heepsing_comp.py ${BEAM_ENERGY} ${MAXEVENTS} ${DATA_Suffix} ${DUMMY_Suffix} ${SIMC_Suffix} ${DATA_RUN_LIST} ${DUMMY_RUN_LIST} ${CSV_FILE} ${Spec}
                    else
                        echo "Skipping python HeeP plotting script step"
                    fi
                elif [ ! -f  "${UTILPATH}/OUTPUT/Analysis/HeeP/${BEAM_ENERGY}_${MAXEVENTS}_HeePSing_${Spec}_Output_Data.root" ]; then
                       python3 ${UTILPATH}/scripts/heep/src/plot_heepsing_comp.py ${BEAM_ENERGY} ${MAXEVENTS} ${DATA_Suffix} ${DUMMY_Suffix} ${SIMC_Suffix} ${DATA_RUN_LIST} ${DUMMY_RUN_LIST} ${CSV_FILE} ${Spec}
                else echo "HeeP coin output plots file already found in ${UTILPATH}/OUTPUT/Analysis/HeeP/ - Skipped python plotting script step"
                fi
                #evince "${UTILPATH}/OUTPUT/Analysis/HeeP/${BEAM_ENERGY}_${MAXEVENTS}_heep_Proton_Analysis_Distributions.pdf" &
	        )
                break;;
            [Nn]* )
                exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

sleep 3

#################################################################################################################################################

# Section for HeeP physics ploting script for beam energy
elif [[ $b_flag = "true" ]]; then
    while true; do
        read -p "Do you wish to analyse and plot the ROOT files for beam energy ${BEAM_ENERGY}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                i=-1
                (
                ##Reads in input file##
                while IFS='' read -r line || [[ -n "$line" ]]; do
                    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                    echo "Run number read from file: $line"
                    echo ""
                    cd "${UTILPATH}/scripts/heep/src/"
                    python3 singyield_heep.py ${ROOTPREFIX_CUT} $line ${MAXEVENTS} ${Spec}
                done < "$inputFile"
                cd $REPLAYPATH/OUTPUT/Analysis/HeeP/
                dir_name="runbyrun"
                # Check if the directory exists
                if [ ! -d "$dir_name" ]; then
                   # If it doesn't exist, create it
                   mkdir "$dir_name"
                   echo "Directory '$dir_name' created."
                else
                   echo "Directory '$dir_name' already exists."
                fi
		if [ -e "${BEAM_ENERGY}_${MAXEVENTS}_HeePSing_${Spec}_Analysed_Data.root" ]; then
                     echo "Deleting ${BEAM_ENERGY}_${MAXEVENTS}_HeePSing_${Spec}_Analysed_Data.root"
                     rm ${BEAM_ENERGY}_${MAXEVENTS}_HeePSing_${Spec}_Analysed_Data.root
                fi
                hadd ${BEAM_ENERGY}_${MAXEVENTS}_HeePSing_${Spec}_Analysed_Data.root 1*
                mv 1* $dir_name
                cd $REPLAYPATH
		##Reads in input file for dummy##
                while IFS='' read -r line || [[ -n "$line" ]]; do
                    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                    echo "Run number read from file: $line"
                    echo ""
                    cd "${UTILPATH}/scripts/heep/src/"
                    python3 singyield_heep.py ${ROOTPREFIX_CUT} $line ${MAXEVENTS} ${Spec}
                done < "${inputFile}_Dummy"
		cd $REPLAYPATH/OUTPUT/Analysis/HeeP/
                if [ -e "${BEAM_ENERGY}_${MAXEVENTS}_HeePSing_${Spec}_Analysed_Dummy_Data.root" ]; then
                        echo "Deleting ${BEAM_ENERGY}_${MAXEVENTS}_HeePSing_${Spec}_Analysed_Dummy_Data.root"
                        rm ${BEAM_ENERGY}_${MAXEVENTS}_HeePSing_${Spec}_Analysed_Dummy_Data.root
                fi
                hadd ${BEAM_ENERGY}_${MAXEVENTS}_HeePSing_${Spec}_Analysed_Dummy_Data.root 1*
                mv 1* $dir_name
		cd $REPLAYPATH
                # Section for HeeP physics ploting script
                if [ -f "${UTILPATH}/OUTPUT/Analysis/HeeP/${BEAM_ENERGY}_${MAXEVENTS}_HeePSing_${Spec}_Output_Data.root" ]; then
                    read -p "HeeP coin output plots file already exists, do you want to reprocess it? <Y/N> " option2
                    if [[ $option2 == "y" || $option2 == "Y" || $option2 == "yes" || $option2 == "Yes" ]]; then
                        rm "${UTILPATH}/OUTPUT/Analysis/HeeP/${BEAM_ENERGY}_${MAXEVENTS}_HeePSing_${Spec}_Output_Data.root"
                        echo "Reprocessing"
                        python3 ${UTILPATH}/scripts/heep/src/plot_heepsing_comp.py ${BEAM_ENERGY} ${MAXEVENTS} ${DATA_Suffix} ${DUMMY_Suffix} ${SIMC_Suffix} ${DATA_RUN_LIST} ${DUMMY_RUN_LIST} ${CSV_FILE} ${Spec}
                    else
                        echo "Skipping python HeeP plotting script step"
                    fi
                elif [ ! -f  "${UTILPATH}/OUTPUT/Analysis/HeeP/${BEAM_ENERGY}_${MAXEVENTS}_HeePSing_${Spec}_Output_Data.root" ]; then
                       python3 ${UTILPATH}/scripts/heep/src/plot_heepsing_comp.py ${BEAM_ENERGY} ${MAXEVENTS} ${DATA_Suffix} ${DUMMY_Suffix} ${SIMC_Suffix} ${DATA_RUN_LIST} ${DUMMY_RUN_LIST} ${CSV_FILE} ${Spec}
                else echo "HeeP coin output plots file already found in ${UTILPATH}/OUTPUT/Analysis/HeeP/ - Skipped python plotting script step"
                fi
                #evince "${UTILPATH}/OUTPUT/Analysis/HeeP/${BEAM_ENERGY}_${MAXEVENTS}_heep_Proton_Analysis_Distributions.pdf" &	
	        )
                break;;
            [Nn]* )
                exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done
fi

exit 0
