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

while getopts 'hcmdetpsra' flag; do
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
        echo "    -a, run script to calculate average kinematics and yields for LTSep analysis "
        echo "        comp -> Flag=-a RunType=arg1-(Prod) RunList=arg2 MaxEvents=arg3-(Optional)"
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
        a) a_flag='true' ;;
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

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

cd $REPLAYPATH
# Input Arguments for Cut Implementation Script
inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/PionLT_2021_2022/${RunList}"
PHY_SETTING=$(echo "${RunList}" | awk -F'_' '{print $1 "_" $2 "_" $3 "_" $4 "_" $5}')
ROOTPREFIX_CUT=${ANATYPE}LT_ProdCoin_replay_production

# Input Arguments for Background Subtraction and MMP Cut Determination Script
SIMC_SETTING=$(echo "${RunList}" | awk -F'_' '{print $1 $2 $3 "_" $4 $5}')
DATA_Suffix=ProdCoin_Analysed_Data
DUMMY_Suffix=ProdCoin_Analysed_Dummy_Data
SIMC_Suffix="Prod_Coin_${SIMC_SETTING}"
DATA_RUN_LIST=${PHY_SETTING}
DUMMY_RUN_LIST=${PHY_SETTING}_dummy
#CSV_FILE=PionLT_coin_production_Prod_efficiency_data_2025_10_23
CSV_FILE=PionLT_coin_production_Prod_efficiency_data_2025_11_21

# Input Arguments for t-resolution and t-binning Scripts
PHY_SETTING_C=$(echo "${RunList}" | awk -F'_' '{print $1 "_" $2 "_" $3}')
SIMC_SETTING_C=$(echo "${RunList}" | awk -F'_' '{print $1 $2 $3}')
SIMC_Suffix_C="Prod_Coin_${SIMC_SETTING_C}"
RUN_LIST_C=${PHY_SETTING_C}_Runlist

# Input Arguments for avergae kinematics and yields calculation Script
DATA_YIELD_CSV=Physics_Data_Yield
DATA_AVG_KIN_CSV=Physics_Avg_Data_Kinematics
SIMC_YIELD_CSV=Physics_SIMC_Yield
SIMC_AVG_KIN_CSV=Physics_Avg_SIMC_Kinematics

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Create directories to store CSVs if it doesn't exist
AVG_KIN_DIR=${UTILPATH}/LTSep_CSVs/avg_kinematics_csv/${PHY_SETTING_C}_std
RATIO_DIR=${UTILPATH}/LTSep_CSVs/datasimc_ratios_csv/${PHY_SETTING_C}_std
DIAMOND_DIR=${UTILPATH}/LTSep_CSVs/diamond_cut_csv/${PHY_SETTING_C}_std
INPUT_DIR=${UTILPATH}/LTSep_CSVs/ltsep_input_csv/${PHY_SETTING_C}_std
MM_OFFSET_DIR=${UTILPATH}/LTSep_CSVs/mm_offset_cut_csv/${PHY_SETTING_C}_std
PHYSICS_YIELDS_DIR=${UTILPATH}/LTSep_CSVs/physics_yields_csv/${PHY_SETTING_C}_std
SIMC_YIELDS_DIR=${UTILPATH}/LTSep_CSVs/simc_yields_csv/${PHY_SETTING_C}_std
T_BINNING_DIR=${UTILPATH}/LTSep_CSVs/t_binning_csv/${PHY_SETTING_C}_std
T_RESOLUTION_DIR=${UTILPATH}/LTSep_CSVs/t_resolution_csv/${PHY_SETTING_C}_std

directories=("${AVG_KIN_DIR}" "${RATIO_DIR}" "${DIAMOND_DIR}" "${INPUT_DIR}" "${MM_OFFSET_DIR}" "${PHYSICS_YIELDS_DIR}" "${SIMC_YIELDS_DIR}" "${T_BINNING_DIR}" "${T_RESOLUTION_DIR}")
# Create all directories
for dir in "${directories[@]}"; do
    if [ ! -d "$dir" ]; then
        echo "Creating directory: $dir"
        mkdir -p "$dir"
    fi
done

# Moving SIMC ROOT files to required standard directory
SIMC_Suffix_File=Prod_Coin_${SIMC_SETTING_C}
# Define source and destination directories
simc_dirs=("${VOLATILEPATH}/worksim" "${VOLATILEPATH}/OUTPUT/Analysis/SIMC")
# Move SIMC files to standard directories
for base_dir in "${simc_dirs[@]}"; do
    dest_dir="${base_dir}/${PHY_SETTING_C}_std"
    # Create destination directory if it doesn't exist
    if [ ! -d "$dest_dir" ]; then
        mkdir -p "$dest_dir"
    fi
    # Only move if destination does NOT already have files
    if ! ls ${dest_dir}/${SIMC_Suffix_File}* 1> /dev/null 2>&1; then
        if ls ${base_dir}/${SIMC_Suffix_File}* 1> /dev/null 2>&1; then
            mv ${base_dir}/${SIMC_Suffix_File}* "$dest_dir/"
        fi
    fi
done

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
                dir_name="${PHY_SETTING_C}_runbyrun"
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
                dcut_dir_name="diamondcut_check"
                cd $REPLAYPATH/OUTPUT/Analysis/PionLT/
                if [ ! -d "$dcut_dir_name" ]; then
                    mkdir "$dcut_dir_name"
                    echo "Directory '$dcut_dir_name' created."
                else
                    echo "Directory '$dcut_dir_name' already exists."
                fi

                mv *Diamondcut* "$dcut_dir_name"
#                cd "$REPLAYPATH/OUTPUT/Analysis/PionLT/$mm_dir_name"
                echo "Output files moved to $dcut_dir_name"
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
        read -p "Do you wish to do t-resolution for physics setting ${PHY_SETTING_C}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                # Processing logic for Pion physics plotting
                output_file="${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING_C}_${MAXEVENTS}_ProdCoin_tresolution_Data.root"

                if [ -f "$output_file" ]; then
                    read -p "Pion coin output plots file already exists, do you want to reprocess it? <Y/N> " option2
                    if [[ $option2 =~ ^[Yy](es)?$ ]]; then
                        rm "$output_file"
                        echo "Reprocessing"
                        python3 ${UTILPATH}/scripts/yield/src/pion_simc_t_resol.py ${PHY_SETTING_C} ${SIMC_SETTING_C}
                    else
                        echo "Skipping python t-resolution script step"
                    fi
                else
                    python3 ${UTILPATH}/scripts/yield/src/pion_simc_t_resol.py ${PHY_SETTING_C} ${SIMC_SETTING_C}
                fi

                # Directory setup
                tresl_dir_name="t_binning"
                cd $REPLAYPATH/OUTPUT/Analysis/PionLT/
                if [ ! -d "$tresl_dir_name" ]; then
                    mkdir "$tresl_dir_name"
                    echo "Directory '$tresl_dir_name' created."
                else
                    echo "Directory '$tresl_dir_name' already exists."
                fi

                mv *tresolution* "$tresl_dir_name"
#                cd "$REPLAYPATH/OUTPUT/Analysis/PionLT/$mm_dir_name"
                echo "Output files moved to $tresl_dir_name"
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
        read -p "Do you wish to do t-binning for physics setting ${PHY_SETTING_C}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                # Processing logic for Pion physics plotting
                output_file="${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING_C}_${MAXEVENTS}_ProdCoin_tbinning_Data.root"

                if [ -f "$output_file" ]; then
                    read -p "Pion coin output plots file already exists, do you want to reprocess it? <Y/N> " option2
                    if [[ $option2 =~ ^[Yy](es)?$ ]]; then
                        rm "$output_file"
                        echo "Reprocessing"
                        python3 ${UTILPATH}/scripts/yield/src/pion_yield_tbinning.py ${PHY_SETTING_C} ${MAXEVENTS} ${DATA_Suffix} ${DUMMY_Suffix} ${RUN_LIST_C}
                    else
                        echo "Skipping python t-binning script step"
                    fi
                else
                    python3 ${UTILPATH}/scripts/yield/src/pion_yield_tbinning.py ${PHY_SETTING_C} ${MAXEVENTS} ${DATA_Suffix} ${DUMMY_Suffix} ${RUN_LIST_C}
                fi

                # Directory setup
                t_dir_name="t_binning"
                cd $REPLAYPATH/OUTPUT/Analysis/PionLT/
                if [ ! -d "$t_dir_name" ]; then
                    mkdir "$t_dir_name"
                    echo "Directory '$t_dir_name' created."
                else
                    echo "Directory '$t_dir_name' already exists."
                fi

                mv *tbinning* "$t_dir_name"
#                cd "$REPLAYPATH/OUTPUT/Analysis/PionLT/$mm_dir_name"
                echo "Output files moved to $t_dir_name"
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
        read -p "Do you wish to do phi binning and calculate yields for physics setting ${PHY_SETTING}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                i=-1
                (
                # Section for Pion physics ploting script
                if [ -f "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING}_${MAXEVENTS}_ProdCoin_Yield_Data.root" ]; then
                    read -p "Pion coin output plots file already exists, do you want to reprocess it? <Y/N> " option2
                    if [[ $option2 == "y" || $option2 == "Y" || $option2 == "yes" || $option2 == "Yes" ]]; then
                        rm "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING}_${MAXEVENTS}_ProdCoin_Yield_Data.root"
                        echo "Reprocessing"
                        python3 ${UTILPATH}/scripts/yield/src/pion_physics_yield.py ${PHY_SETTING} ${MAXEVENTS} ${DATA_Suffix} ${DUMMY_Suffix} ${DATA_RUN_LIST} ${DUMMY_RUN_LIST} ${CSV_FILE}
                    else
                        echo "Skipping phi-binning and yield calculation script step"
                    fi
                elif [ ! -f  "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING}_${MAXEVENTS}_ProdCoin_Yield_Data.root" ]; then
                       python3 ${UTILPATH}/scripts/yield/src/pion_physics_yield.py ${PHY_SETTING} ${MAXEVENTS} ${DATA_Suffix} ${DUMMY_Suffix} ${DATA_RUN_LIST} ${DUMMY_RUN_LIST} ${CSV_FILE}
                else echo "Pion coin output plots file already found in ${UTILPATH}/OUTPUT/Analysis/PionLT/ - Skipped python plotting script step"
                fi
	        )
                # Directory setup
                phy_dir_name="data_yield"
                cd $REPLAYPATH/OUTPUT/Analysis/PionLT/
                if [ ! -d "$phy_dir_name" ]; then
                    mkdir "$phy_dir_name"
                    echo "Directory '$phy_dir_name' created."
                else
                    echo "Directory '$phy_dir_name' already exists."
                fi

                mv *Yield* "$phy_dir_name"
#                cd "$REPLAYPATH/OUTPUT/Analysis/PionLT/$mm_dir_name"
                echo "Output files moved to $phy_dir_name"
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
        read -p "Do you wish to do t & phi binning and calculate yields for simc setting ${PHY_SETTING}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                i=-1
                (
                # Section for Pion simc ploting script
                if [ -f "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING}_${MAXEVENTS}_ProdCoin_Yield_Data.root" ]; then
                    read -p "Pion coin output plots file already exists, do you want to reprocess it? <Y/N> " option2
                    if [[ $option2 == "y" || $option2 == "Y" || $option2 == "yes" || $option2 == "Yes" ]]; then
                        echo "Reprocessing"
                        python3 ${UTILPATH}/scripts/yield/src/pion_simc_yield.py ${PHY_SETTING} ${MAXEVENTS} ${SIMC_Suffix}
                    else
                        echo "Skipping phi-binning and yield calculation script step"
                    fi
                elif [ ! -f  "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING}_${MAXEVENTS}_ProdCoin_Yield_Data.root" ]; then
                       python3 ${UTILPATH}/scripts/yield/src/pion_simc_yield.py ${PHY_SETTING} ${MAXEVENTS} ${SIMC_Suffix}
                else echo "Pion coin output plots file already found in ${UTILPATH}/OUTPUT/Analysis/PionLT/ - Skipped python plotting script step"
                fi
	        )
                # Directory setup for ROOT files
                cd $REPLAYPATH/OUTPUT/Analysis/PionLT/
                simc_dir_yield="simc_yield"
                simc_dir_iter="${PHY_SETTING_C}_std"
                if [ ! -d "$simc_dir_yield/$simc_dir_iter" ]; then
                    mkdir -p "$simc_dir_yield/$simc_dir_iter"
                    echo "Directory '$simc_dir_yield/$simc_dir_iter' created."
                else
                    echo "Directory '$simc_dir_yield/$simc_dir_iter' already exists."
                fi
                mv *Yield* "$simc_dir_yield/$simc_dir_iter"
                echo "Output files moved to $simc_dir_yield/$simc_dir_iter"
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
        read -p "Do you wish to do Data & SIMC comparison and plotting for physics setting ${PHY_SETTING}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                i=-1
                (
                # Section for Pion physics ploting script
                if [ -f "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING}_${MAXEVENTS}_ProdCoin_Yield_Data.root" ]; then
                    read -p "Pion coin output plots file already exists, do you want to reprocess it? <Y/N> " option2
                    if [[ $option2 == "y" || $option2 == "Y" || $option2 == "yes" || $option2 == "Yes" ]]; then
                        rm "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING}_${MAXEVENTS}_ProdCoin_Yield_Data.root"
                        echo "Reprocessing"
                        python3 ${UTILPATH}/scripts/yield/src/pion_physics_ratio_comp.py ${PHY_SETTING} ${MAXEVENTS} ${DATA_Suffix} ${DUMMY_Suffix} ${SIMC_Suffix} ${DATA_RUN_LIST} ${DUMMY_RUN_LIST} ${CSV_FILE}
                    else
                        echo "Skipping  Data & SIMC comparison and plotting script step"
                    fi
                elif [ ! -f  "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING}_${MAXEVENTS}_ProdCoin_Yield_Data.root" ]; then
                       python3 ${UTILPATH}/scripts/yield/src/pion_physics_ratio_comp.py ${PHY_SETTING} ${MAXEVENTS} ${DATA_Suffix} ${DUMMY_Suffix} ${SIMC_Suffix} ${DATA_RUN_LIST} ${DUMMY_RUN_LIST} ${CSV_FILE}
                else echo "Pion coin output plots file already found in ${UTILPATH}/OUTPUT/Analysis/PionLT/ - Skipped python plotting script step"
                fi
	        )
                # Directory setup for ROOT files
                cd $REPLAYPATH/OUTPUT/Analysis/PionLT/
                ratio_dir="ratios"
                ratio_dir_iter="${PHY_SETTING_C}_std"
                if [ ! -d "$ratio_dir/$ratio_dir_iter" ]; then
                    mkdir -p "$ratio_dir/$ratio_dir_iter"
                    echo "Directory '$ratio_dir/$ratio_dir_iter' created."
                else
                    echo "Directory '$ratio_dir/$ratio_dir_iter' already exists."
                fi
                mv *Ratio* "$ratio_dir/$ratio_dir_iter"
                echo "Output files moved to $ratio_dir/$ratio_dir_iter"
                break;;
            [Nn]* )
                exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

sleep 3

#################################################################################################################################
# Section for Average Kinematics and Yield Calculation Script
elif [[ $a_flag == "true" ]]; then
    while true; do
        read -p "Do you wish to calculate average kinematics and yields for physics setting ${PHY_SETTING_C}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                i=-1
                (
                # Section for Pion physics ploting script
                if [ -f "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING_C}_${MAXEVENTS}_ProdCoin_Yield_Data.root" ]; then
                    read -p "Pion coin output plots file already exists, do you want to reprocess it? <Y/N> " option2
                    if [[ $option2 == "y" || $option2 == "Y" || $option2 == "yes" || $option2 == "Yes" ]]; then
                        rm "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING_C}_${MAXEVENTS}_ProdCoin_Yield_Data.root"
                        echo "Reprocessing"
                        python3 ${UTILPATH}/scripts/yield/src/pion_yield_avg_kin.py ${PHY_SETTING_C} ${MAXEVENTS} ${DATA_YIELD_CSV} ${DATA_AVG_KIN_CSV} ${SIMC_YIELD_CSV} ${SIMC_AVG_KIN_CSV}
                    else
                        echo "Skipping  Data & SIMC comparison and plotting script step"
                    fi
                elif [ ! -f  "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING_C}_${MAXEVENTS}_ProdCoin_Yield_Data.root" ]; then
                       python3 ${UTILPATH}/scripts/yield/src/pion_yield_avg_kin.py ${PHY_SETTING_C} ${MAXEVENTS} ${DATA_YIELD_CSV} ${DATA_AVG_KIN_CSV} ${SIMC_YIELD_CSV} ${SIMC_AVG_KIN_CSV}
                else echo "Pion coin output plots file already found in ${UTILPATH}/OUTPUT/Analysis/PionLT/ - Skipped python plotting script step"
                fi
	        )
                # Directory setup for ROOT files
                cd $REPLAYPATH/OUTPUT/Analysis/PionLT/
                avgkin_dir="avg_kin"
                avg_kin_dir_iter="${PHY_SETTING_C}_std"
                if [ ! -d "$avgkin_dir/$avg_kin_dir_iter" ]; then
                    mkdir -p "$avgkin_dir/$avg_kin_dir_iter"
                    echo "Directory '$avgkin_dir/$avg_kin_dir_iter' created."
                else
                    echo "Directory '$avgkin_dir/$avg_kin_dir_iter' already exists."
                fi
                mv *avgkin* "$avgkin_dir/$avg_kin_dir_iter"
                echo "Output files moved to $avgkin_dir/$avg_kin_dir_iter"
                break;;
            [Nn]* )
                exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done
fi

exit 0