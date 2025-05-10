#! /bin/bash

#
# Description:
# ================================================================
# Time-stamp: "2025-03-11 12:12:00 junaid"
# ================================================================
#
# Author:  Muhammad Junaid <mjo147@uregina.ca>
#
# Copyright (c) junaid
#
##################################################################################
# Created - 10/July/2021, Author - Muhammad Junaid, University of Regina, Canada
##################################################################################

while getopts 'cmdetpsr' flag; do
    case "${flag}" in
        h)
        echo "-------------------------------------------------------------------"
        echo "./pion_prodYield_kinematic_analysis.sh -{flags} {variable arguments, see help}"
        echo "-------------------------------------------------------------------"
        echo
        echo "The following flags can be called for the physics analysis..."         
        echo "    -h, help"
        echo "    -c, run cut analyser script to trim down the root file for further analysis"
        echo "        cut -> Flag=-c RunType=arg1-(Prod) RunList=arg2 MaxEvents=arg3-(Optional)"        
        echo "    -m, run missing mass check script to determine the cut for missing mass"
        echo "        mmcut -> Flag=-m RunType=arg1-(Prod) RunList=arg2 MaxEvents=arg3-(Optional)"
        echo "    -d, run diamond cut check script to determine the diamond cut for data"
        echo "        dcut -> Flag=-d RunType=arg1-(Prod) RunList=arg2 MaxEvents=arg3-(Optional)"
        echo "    -e, run t-binning script to determine the t-resolution for data"
        echo "        tresl -> Flag=-r RunType=arg1-(Prod) Setting=arg2 MaxEvents=arg3-(Optional)"
        echo "    -t, run t-binning script to determine the t-binning for data"
        echo "        tbin -> Flag=-t RunType=arg1-(Prod) Setting=arg2 MaxEvents=arg3-(Optional)"
        echo "    -p, run phi-binning and physics yield script to determine the physics yields"
        echo "        phibin_pyield -> Flag=-y RunType=arg1-(Prod) RunList=arg2 MaxEvents=arg3-(Optional)"     
        echo "    -s, run simc yield script to determine simc yields"
        echo "        simc -> Flag=-s RunType=arg1-(Prod) RunList=arg2 MaxEvents=arg3-(Optional)"
        echo "    -r, run plotting script to compare Data/SIMC and determine ratios"
        echo "        comp -> Flag=-p RunType=arg1-(Prod) RunList=arg2 MaxEvents=arg3-(Optional)"
        exit 0
        ;;
        c) c_flag='true' ;;
        m) m_flag='true' ;;
        d) d_flag='true' ;;
        e) e_flag='true' ;;
        t) t_flag='true' ;;
        p) p_flag='true' ;;
        s) s_flag='true' ;;
        r) r_flag='true' ;;
        *) print_usage
        exit 1 ;;
    esac
done

echo "Starting analysis of Pion events"
echo "I take as arguments the -flag, physics setting, runlist and number of events!"

# Input params - beam energy, run list and max number of events
RunType=$2
if [[ -z "$2" ]]; then
    echo "No RunType provided. Please provide type of Analysis 'ProdCoin or Prod'!"
fi

RunList=$3
if [[ -z "$3" ]]; then
    echo "I need a runlist"
    echo "Please provide a run list as input"
fi

MAXEVENTS=$4
if [[ -z "$4" ]]; then
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
# Input Arguments for Cut Implementation Script
inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/PionLT_2021_2022/${RunList}"
PHY_SETTING=$(echo "${RunList}" | awk -F'_' '{print $1 "_" $2 "_" $3 "_" $4 "_" $5}')
ROOTPREFIX_CUT=${ANATYPE}LT_ProdCoin_replay_production

# Input Arguments for Background Subtraction and MMP Cut Determination Script
DATA_Suffix=ProdCoin_Analysed_Data
DUMMY_Suffix=ProdCoin_Analysed_Dummy_Data
SIMC_Suffix="Prod_Coin_$(echo "${RunList}" | awk -F'_' '{print $1 $2 $3 "_" $4 $5}')"
DATA_RUN_LIST=${PHY_SETTING}
DUMMY_RUN_LIST=${PHY_SETTING}_dummy
CSV_FILE=PionLT_coin_production_Prod_efficiency_data_2025_03_08

# Input Arguments for t-resolution and t-binning Scripts
PHY_SETTING_tbin=$(echo "${RunList}" | awk -F'_' '{print $1 "_" $2 "_" $3}')
DATA_Suffix_tbin=$(echo "${RunList}" | awk -F'_' '{print $1 "_" $2 "_" $3}')
SIMC_Suffix_tbin=$(echo "${RunList}" | awk -F'_' '{print $1 $2 $3}')
RUN_LIST_tbin=${PHY_SETTING_tbin}_Runlist

################################################################################################################################                                                                                   

# Section for Pion Physics Analysis Script
if [[ $c_flag == "true" ]]; then
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
                    cd "${UTILPATH}/scripts/yield/src/"
                    python3 pion_yield_cuts.py ${ROOTPREFIX_CUT} $line ${MAXEVENTS}
                done < "$inputFile"
                cd $REPLAYPATH/OUTPUT/Analysis/PionLT/
                dir_name="runbyrun"
                # Check if the directory exists
                if [ ! -d "$dir_name" ]; then
                   # If it doesn't exist, create it
                   mkdir "$dir_name"
                   echo "Directory '$dir_name' created."
                else
                   echo "Directory '$dir_name' already exists."
                fi
                if [[ "${RunList,,}" == *"dummy"* ]]; then
                    if [ -e "${PHY_SETTING}_${MAXEVENTS}_ProdCoin_Analysed_Dummy_Data.root" ]; then
                        echo "Deleting ${PHY_SETTING}_${MAXEVENTS}_ProdCoin_Analysed_Dummy_Data.root"
    		        rm ${PHY_SETTING}_${MAXEVENTS}_ProdCoin_Analysed_Dummy_Data.root
                    fi  
		    hadd ${PHY_SETTING}_${MAXEVENTS}_ProdCoin_Analysed_Dummy_Data.root 1*
		else
                    if [ -e "${PHY_SETTING}_${MAXEVENTS}_ProdCoin_Analysed_Data.root" ]; then
			echo "Deleting ${PHY_SETTING}_${MAXEVENTS}_ProdCoin_Analysed_Data.root"
			rm ${PHY_SETTING}_${MAXEVENTS}_ProdCoin_Analysed_Data.root
                    fi
		    hadd ${PHY_SETTING}_${MAXEVENTS}_ProdCoin_Analysed_Data.root 1*
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
# Section for MM Cut Check script
elif [[ $m_flag == "true" ]]; then
    while true; do
        read -p "Do you wish to do Missing Mass cut check for physics setting ${PHY_SETTING}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                # Processing logic for Pion physics plotting
                output_file="${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING}_${MAXEVENTS}_ProdCoin_MMcut_Data.root"
                
                if [ -f "$output_file" ]; then
                    read -p "Pion coin output plots file already exists, do you want to reprocess it? <Y/N> " option2
                    if [[ $option2 =~ ^[Yy](es)?$ ]]; then
                        rm "$output_file"
                        echo "Reprocessing"
                        python3 ${UTILPATH}/scripts/yield/src/pion_yield_MMcut.py ${PHY_SETTING} ${MAXEVENTS} ${DATA_Suffix} ${DUMMY_Suffix} ${SIMC_Suffix} ${DATA_RUN_LIST} ${DUMMY_RUN_LIST} ${CSV_FILE}
                    else
                        echo "Skipping python Missing Mass cut script step"
                    fi
                else
                    python3 ${UTILPATH}/scripts/yield/src/pion_yield_MMcut.py ${PHY_SETTING} ${MAXEVENTS} ${DATA_Suffix} ${DUMMY_Suffix} ${SIMC_Suffix} ${DATA_RUN_LIST} ${DUMMY_RUN_LIST} ${CSV_FILE}
                fi

                # Directory setup
                mm_dir_name="mmcut_check"
                cd $REPLAYPATH/OUTPUT/Analysis/PionLT/
                if [ ! -d "$mm_dir_name" ]; then
                    mkdir "$mm_dir_name"
                    echo "Directory '$mm_dir_name' created."
                else
                    echo "Directory '$mm_dir_name' already exists."
                fi
                
                mv *MMcut* "$mm_dir_name"
#                cd "$REPLAYPATH/OUTPUT/Analysis/PionLT/$mm_dir_name"
                echo "Output files moved to $mm_dir_name"
                break;;
            [Nn]* )
                exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

sleep 3

#################################################################################################################################

# Section for Diamond Cut script
elif [[ $d_flag == "true" ]]; then
    while true; do
        read -p "Do you wish to do Diamond cut check for physics setting ${PHY_SETTING}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                # Processing logic for Pion physics plotting
                output_file="${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING}_${MAXEVENTS}_ProdCoin_MMcut_Data.root"
                
                if [ -f "$output_file" ]; then
                    read -p "Pion coin output plots file already exists, do you want to reprocess it? <Y/N> " option2
                    if [[ $option2 =~ ^[Yy](es)?$ ]]; then
                        rm "$output_file"
                        echo "Reprocessing"
                        python3 ${UTILPATH}/scripts/yield/src/pion_yield_diamondcut.py ${PHY_SETTING} ${MAXEVENTS} ${DATA_Suffix} ${DUMMY_Suffix} ${SIMC_Suffix} ${DATA_RUN_LIST} ${DUMMY_RUN_LIST} ${CSV_FILE}
                    else
                        echo "Skipping python Diamond cut script step"
                    fi
                else
                    python3 ${UTILPATH}/scripts/yield/src/pion_yield_diamondcut.py ${PHY_SETTING} ${MAXEVENTS} ${DATA_Suffix} ${DUMMY_Suffix} ${SIMC_Suffix} ${DATA_RUN_LIST} ${DUMMY_RUN_LIST} ${CSV_FILE}
                fi

                # Directory setup
                mm_dir_name="diamondcut_check"
                cd $REPLAYPATH/OUTPUT/Analysis/PionLT/
                if [ ! -d "$mm_dir_name" ]; then
                    mkdir "$mm_dir_name"
                    echo "Directory '$mm_dir_name' created."
                else
                    echo "Directory '$mm_dir_name' already exists."
                fi
                
                mv *Diamondcut* "$mm_dir_name"
#                cd "$REPLAYPATH/OUTPUT/Analysis/PionLT/$mm_dir_name"
                echo "Output files moved to $mm_dir_name"
                break;;
            [Nn]* )
                exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

sleep 3

#################################################################################################################################

# Section for  t-resolution script
elif [[ $e_flag == "true" ]]; then
    while true; do
        read -p "Do you wish to do t-resolution for physics setting ${PHY_SETTING2}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                # Processing logic for Pion physics plotting
                output_file="${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING2}_${MAXEVENTS}_ProdCoin_tresolution_Data.root"
                
                if [ -f "$output_file" ]; then
                    read -p "Pion coin output plots file already exists, do you want to reprocess it? <Y/N> " option2
                    if [[ $option2 =~ ^[Yy](es)?$ ]]; then
                        rm "$output_file"
                        echo "Reprocessing"
                        python3 ${UTILPATH}/scripts/yield/src/pion_simc_t_resol.py ${PHY_SETTING_tbin} ${SIMC_Suffix_tbin} 
                    else
                        echo "Skipping python t-resolution script step"
                    fi
                else
                    python3 ${UTILPATH}/scripts/yield/src/pion_simc_t_resol.py ${PHY_SETTING_tbin} ${SIMC_Suffix_tbin}
                fi

                # Directory setup
                mm_dir_name="t_binning"
                cd $REPLAYPATH/OUTPUT/Analysis/PionLT/
                if [ ! -d "$mm_dir_name" ]; then
                    mkdir "$mm_dir_name"
                    echo "Directory '$mm_dir_name' created."
                else
                    echo "Directory '$mm_dir_name' already exists."
                fi
                
                mv *tresolution* "$mm_dir_name"
#                cd "$REPLAYPATH/OUTPUT/Analysis/PionLT/$mm_dir_name"
                echo "Output files moved to $mm_dir_name"
                break;;
            [Nn]* )
                exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

sleep 3

#################################################################################################################################

# Section for t-binning script
elif [[ $t_flag == "true" ]]; then
    while true; do
        read -p "Do you wish to do t-binning for physics setting ${PHY_SETTING2}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                # Processing logic for Pion physics plotting
                output_file="${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING2}_${MAXEVENTS}_ProdCoin_tbinning_Data.root"
                
                if [ -f "$output_file" ]; then
                    read -p "Pion coin output plots file already exists, do you want to reprocess it? <Y/N> " option2
                    if [[ $option2 =~ ^[Yy](es)?$ ]]; then
                        rm "$output_file"
                        echo "Reprocessing"
                        python3 ${UTILPATH}/scripts/yield/src/pion_yield_tbinning.py ${PHY_SETTING_tbin} ${MAXEVENTS} ${DATA_Suffix_tbin} ${RUN_LIST_tbin}
                    else
                        echo "Skipping python t-binning script step"
                    fi
                else
                    python3 ${UTILPATH}/scripts/yield/src/pion_yield_tbinning.py ${PHY_SETTING_tbin} ${MAXEVENTS} ${DATA_Suffix_tbin} ${RUN_LIST_tbin}
                fi

                # Directory setup
                mm_dir_name="t_binning"
                cd $REPLAYPATH/OUTPUT/Analysis/PionLT/
                if [ ! -d "$mm_dir_name" ]; then
                    mkdir "$mm_dir_name"
                    echo "Directory '$mm_dir_name' created."
                else
                    echo "Directory '$mm_dir_name' already exists."
                fi
                
                mv *tbinning* "$mm_dir_name"
#                cd "$REPLAYPATH/OUTPUT/Analysis/PionLT/$mm_dir_name"
                echo "Output files moved to $mm_dir_name"
                break;;
            [Nn]* )
                exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

sleep 3

#################################################################################################################################
# Section for phi-bining and physics yield calculation Script
elif [[ $p_flag == "true" ]]; then
    while true; do
        read -p "Do you wish to do phi binning and calculate yields for physics setting ${PHY_SETTING2}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                i=-1
                (
                # Section for Pion physics ploting script
                if [ -f "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING2}_${MAXEVENTS}_ProdCoin_Yield_Data.root" ]; then
                    read -p "Pion coin output plots file already exists, do you want to reprocess it? <Y/N> " option2
                    if [[ $option2 == "y" || $option2 == "Y" || $option2 == "yes" || $option2 == "Yes" ]]; then
                        rm "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING2}_${MAXEVENTS}_ProdCoin_Yield_Data.root"
                        echo "Reprocessing"
                        python3 ${UTILPATH}/scripts/yield/src/pion_physics_yield.py ${PHY_SETTING} ${MAXEVENTS} ${DATA_Suffix} ${DUMMY_Suffix} ${DATA_RUN_LIST} ${DUMMY_RUN_LIST} ${CSV_FILE}
                    else
                        echo "Skipping phi-binning and yield calculation script step"
                    fi
                elif [ ! -f  "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING2}_${MAXEVENTS}_ProdCoin_Yield_Data.root" ]; then
                       python3 ${UTILPATH}/scripts/yield/src/pion_physics_yield.py ${PHY_SETTING} ${MAXEVENTS} ${DATA_Suffix} ${DUMMY_Suffix} ${DATA_RUN_LIST} ${DUMMY_RUN_LIST} ${CSV_FILE}
                else echo "Pion coin output plots file already found in ${UTILPATH}/OUTPUT/Analysis/PionLT/ - Skipped python plotting script step"
                fi
	        )
                # Directory setup
                mm_dir_name="yield"
                cd $REPLAYPATH/OUTPUT/Analysis/PionLT/
                if [ ! -d "$mm_dir_name" ]; then
                    mkdir "$mm_dir_name"
                    echo "Directory '$mm_dir_name' created."
                else
                    echo "Directory '$mm_dir_name' already exists."
                fi
                
                mv *Yield* "$mm_dir_name"
#                cd "$REPLAYPATH/OUTPUT/Analysis/PionLT/$mm_dir_name"
                echo "Output files moved to $mm_dir_name"
                break;;
            [Nn]* )
                exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

sleep 3

#################################################################################################################################
# Section for phi-bining and simc yield calculation Script
elif [[ $s_flag == "true" ]]; then
    while true; do
        read -p "Do you wish to do t & phi binning and calculate yields for simc setting ${PHY_SETTING2}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                i=-1
                (
                # Section for Pion simc ploting script
                if [ -f "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING2}_${MAXEVENTS}_ProdCoin_Yield_Data.root" ]; then
                    read -p "Pion coin output plots file already exists, do you want to reprocess it? <Y/N> " option2
                    if [[ $option2 == "y" || $option2 == "Y" || $option2 == "yes" || $option2 == "Yes" ]]; then
                        rm "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING2}_${MAXEVENTS}_ProdCoin_Yield_Data.root"
                        echo "Reprocessing"
                        python3 ${UTILPATH}/scripts/yield/src/pion_simc_yield.py ${PHY_SETTING} ${MAXEVENTS} ${SIMC_Suffix}
                    else
                        echo "Skipping phi-binning and yield calculation script step"
                    fi
                elif [ ! -f  "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING2}_${MAXEVENTS}_ProdCoin_Yield_Data.root" ]; then
                       python3 ${UTILPATH}/scripts/yield/src/pion_simc_yield.py ${PHY_SETTING} ${MAXEVENTS} ${SIMC_Suffix}
                else echo "Pion coin output plots file already found in ${UTILPATH}/OUTPUT/Analysis/PionLT/ - Skipped python plotting script step"
                fi
	        )
                # Directory setup
                mm_dir_name="yield"
                cd $REPLAYPATH/OUTPUT/Analysis/PionLT/
                if [ ! -d "$mm_dir_name" ]; then
                    mkdir "$mm_dir_name"
                    echo "Directory '$mm_dir_name' created."
                else
                    echo "Directory '$mm_dir_name' already exists."
                fi
                
                mv *Yield* "$mm_dir_name"
#                cd "$REPLAYPATH/OUTPUT/Analysis/PionLT/$mm_dir_name"
                echo "Output files moved to $mm_dir_name"
                break;;
            [Nn]* )
                exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

sleep 3

#################################################################################################################################
# Section for Data and SIMC comparison Script
elif [[ $r_flag == "true" ]]; then
    while true; do
        read -p "Do you wish to do Data & SIMC comparison and plotting for physics setting ${PHY_SETTING2}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                i=-1
                (
                # Section for Pion physics ploting script
                if [ -f "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING2}_${MAXEVENTS}_ProdCoin_Yield_Data.root" ]; then
                    read -p "Pion coin output plots file already exists, do you want to reprocess it? <Y/N> " option2
                    if [[ $option2 == "y" || $option2 == "Y" || $option2 == "yes" || $option2 == "Yes" ]]; then
                        rm "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING2}_${MAXEVENTS}_ProdCoin_Yield_Data.root"
                        echo "Reprocessing"
                        python3 ${UTILPATH}/scripts/yield/src/pion_physics_ratio_comp.py ${PHY_SETTING} ${MAXEVENTS} ${DATA_Suffix} ${DUMMY_Suffix} ${SIMC_Suffix} ${DATA_RUN_LIST} ${DUMMY_RUN_LIST} ${CSV_FILE}
                    else
                        echo "Skipping  Data & SIMC comparison and plotting script step"
                    fi
                elif [ ! -f  "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING2}_${MAXEVENTS}_ProdCoin_Yield_Data.root" ]; then
                       python3 ${UTILPATH}/scripts/yield/src/pion_physics_ratio_comp.py ${PHY_SETTING} ${MAXEVENTS} ${DATA_Suffix} ${DUMMY_Suffix} ${SIMC_Suffix} ${DATA_RUN_LIST} ${DUMMY_RUN_LIST} ${CSV_FILE}
                else echo "Pion coin output plots file already found in ${UTILPATH}/OUTPUT/Analysis/PionLT/ - Skipped python plotting script step"
                fi
	        )
                break;;
            [Nn]* )
                exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

fi

exit 0