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
        echo "        cut -> RunNum=arg1 MaxEvents=arg2"        
        echo "    -p, run plotting script to save output as pdf"
        echo "        plot -> RunNum=arg1 MaxEvents=arg2"  
        echo "    -b, run cut analyser and plotting script to save output as pdf for beamenergy"
        echo "        plot -> RunNum=arg1 MaxEvents=arg2"  
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
echo "I take as arguments the run number and number of events!"

# Input params - run number and max number of events
RunType=$2
if [[ -z "$2" ]]; then
    echo "Only Run Number entered...I'll assume, its HeePCoin run"
    RunType=HeePCoin
fi
RUNNUMBER=$3
if [[ -z "$3" ]]; then
    echo "I need a run number"
    echo "Please provide a run number as input"
fi
MAXEVENTS=$4
if [[ -z "$4" ]]; then
    echo "Only Run Number entered...I'll assume -1 (all) events!" 
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
ROOTPREFIX_CUT=${ANATYPE}LT_HeePCoin_replay_production

#################################################################################################################################################

if [[ $c_flag = "true" ]]; then
    while true; do
        read -p "Do you wish to analyse root files with run number ${RUNNUMBER}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                i=-1
                (
                # Section for HeeP analysis script
                if [ -f "${UTILPATH}/OUTPUT/Analysis/HeeP/${RUNNUMBER}_${MAXEVENTS}_Analysed_Data.root" ]; then
                      read -p "HeeP coin production analysed file already exists, do you want to reprocess it? <Y/N> " option1
                if [[ $option1 == "y" || $option1 == "Y" || $option1 == "yes" || $option1 == "Yes" ]]; then
	              rm "${UTILPATH}/OUTPUT/Analysis/HeeP/${RUNNUMBER}_${MAXEVENTS}_Analysed_Data.root"
	        echo "Reprocessing"
	        python3 ${UTILPATH}/scripts/heep/src/coinyield_heep.py ${ROOTPREFIX_CUT} ${RUNNUMBER} ${MAXEVENTS}
                else
	              echo "Skipping python analysis script step"
                fi
                elif [ ! -f "${UTILPATH}/OUTPUT/Analysis/HeeP/${RUNNUMBER}_${MAXEVENTS}_Analysed_Data.root" ]; then
	              python3 ${UTILPATH}/scripts/heep/src/coinyield_heep.py ${ROOTPREFIX_CUT} ${RUNNUMBER} ${MAXEVENTS}
                else echo "Analysed root file already found in ${UTILPATH}/OUTPUT/Analysis/HeeP/ - Skipped python analyzer script step"
                fi
                )
                break;;
            [Nn]* )
                exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

sleep 3

##################################################################################################################################

elif [[ $p_flag = "true" ]]; then
    while true; do
        read -p "Do you wish to plot the ROOT files with runlist ${RUNNUMBER}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                i=-1
                (
                # Section for HeeP physics ploting script
                if [ -f "${UTILPATH}/OUTPUT/Analysis/HeeP/${RUNNUMBER}_${MAXEVENTS}_Output_Data.root" ]; then
                   read -p "HeeP coin output plots file already exists, do you want to reprocess it? <Y/N> " option2
                    if [[ $option2 == "y" || $option2 == "Y" || $option2 == "yes" || $option2 == "Yes" ]]; then
	                rm "${UTILPATH}/OUTPUT/Analysis/HeeP/${RUNNUMBER}_${MAXEVENTS}_Output_Data.root"
	           echo "Reprocessing"
	           python3 ${UTILPATH}/scripts/heep/src/plot_coin_runnumber.py Analysed_Data ${RUNNUMBER} ${MAXEVENTS} 
                else
	           echo "Skipping python HeeP plotting script step"
                fi
                elif [ ! -f  "${UTILPATH}/OUTPUT/Analysis/HeeP/${RUNNUMBER}_${MAXEVENTS}_Output_Data.root" ]; then
	           python3 ${UTILPATH}/scripts/heep/src/plot_coin_runnumber.py Analysed_Data ${RUNNUMBER} ${MAXEVENTS} 
                else echo "HeeP coin output plots file already found in ${UTILPATH}/OUTPUT/Analysis/HeeP/ - Skipped python plotting script step"
                fi
                #evince "${UTILPATH}/OUTPUT/Analysis/HeeP/${RUNNUMBER}_${MAXEVENTS}_heep_Proton_Analysis_Distributions.pdf" &
                )
                break;;
            [Nn]* )
                exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

sleep 3

################################################################################################################################

elif [[ $b_flag = "true" ]]; then
    while true; do
        read -p "Do you wish to analyse and plot the ROOT files for beam energy ${BEAM_ENERGY}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                i=-1
                (
                # Section for HeeP analysis script
                if [ -f "${UTILPATH}/OUTPUT/Analysis/HeeP/${RUNNUMBER}_${MAXEVENTS}_Analysed_Data.root" ]; then
                      read -p "HeeP coin production analysed file already exists, do you want to reprocess it? <Y/N> " option1
                if [[ $option1 == "y" || $option1 == "Y" || $option1 == "yes" || $option1 == "Yes" ]]; then
                      rm "${UTILPATH}/OUTPUT/Analysis/HeeP/${RUNNUMBER}_${MAXEVENTS}_Analysed_Data.root"
                echo "Reprocessing"
                python3 ${UTILPATH}/scripts/heep/src/coinyield_heep.py ${ROOTPREFIX_CUT} ${RUNNUMBER} ${MAXEVENTS}
                else
                      echo "Skipping python analysis script step"
                fi
                elif [ ! -f "${UTILPATH}/OUTPUT/Analysis/HeeP/${RUNNUMBER}_${MAXEVENTS}_Analysed_Data.root" ]; then
                      python3 ${UTILPATH}/scripts/heep/src/coinyield_heep.py ${ROOTPREFIX_CUT} ${RUNNUMBER} ${MAXEVENTS}
                else echo "Analysed root file already found in ${UTILPATH}/OUTPUT/Analysis/HeeP/ - Skipped python analyzer script step"
                fi
                # Section for HeeP physics ploting script
                if [ -f "${UTILPATH}/OUTPUT/Analysis/HeeP/${RUNNUMBER}_${MAXEVENTS}_Output_Data.root" ]; then
                   read -p "HeeP coin output plots file already exists, do you want to reprocess it? <Y/N> " option2
                    if [[ $option2 == "y" || $option2 == "Y" || $option2 == "yes" || $option2 == "Yes" ]]; then
                        rm "${UTILPATH}/OUTPUT/Analysis/HeeP/${RUNNUMBER}_${MAXEVENTS}_Output_Data.root"
                   echo "Reprocessing"
                   python3 ${UTILPATH}/scripts/heep/src/plot_coin_runnumber.py Analysed_Data ${RUNNUMBER} ${MAXEVENTS}
                else
                   echo "Skipping python HeeP plotting script step"
                fi
                elif [ ! -f  "${UTILPATH}/OUTPUT/Analysis/HeeP/${RUNNUMBER}_${MAXEVENTS}_Output_Data.root" ]; then
                   python3 ${UTILPATH}/scripts/heep/src/plot_coin_runnumber.py Analysed_Data ${RUNNUMBER} ${MAXEVENTS}
                else echo "HeeP coin output plots file already found in ${UTILPATH}/OUTPUT/Analysis/HeeP/ - Skipped python plotting script step"
                fi
                #evince "${UTILPATH}/OUTPUT/Analysis/HeeP/${RUNNUMBER}_${MAXEVENTS}_heep_Proton_Analysis_Distributions.pdf" &
         	)
                break;;
            [Nn]* )
                exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

fi
exit 0
