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
        echo "        cut -> RunList=arg1 Beam Energy=arg2 MaxEvents=arg3"        
        echo "    -p, run plotting script to save output as pdf"
        echo "        plot -> RunList=arg1 Beam Energy=arg2 MaxEvents=arg3"  
        echo "    -b, run plotting script to save output as pdf for beamenergy"
        echo "        plot -> RunList=arg1 Beam Energy=arg2 MaxEvents=arg3"  
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
BEAMENERGY=$3
if [[ -z "$3" ]]; then
    echo "Only Run list entered...I'll assume no beam energy value required for cut analysis" 
    BEAMENERGY=0
fi
RunList=$2
if [[ -z "$2" ]]; then
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
    PATHFILE_INFO=`python3 /u/home/${USER}/.local/lib/python3.4/site-packages/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
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

################################################################################################################################                                                                                   

# Section for HeeP analysis script
ROOTPREFIX_CUT=${ANATYPE}LT_HeePCoin_replay_production

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
                    python3 coinyield.py ${ROOTPREFIX_CUT} $line ${MAXEVENTS}
                done < "$inputFile"
                )
                break;;
            [Nn]* )
                exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

##################################################################################################################################

# Section for HeeP physics ploting script
#ROOTPREFIX_PLOT=Analysed_Data

elif [[ $p_flag = "true" ]]; then
    while true; do
        read -p "Do you wish to plot the ROOT files with runlist ${inputFile}? (Please answer yes or no) " yn
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
                    python3 plot_coin.py Analysed_Data $line ${MAXEVENTS}
                done < "$inputFile"
                )
                break;;
            [Nn]* )
                exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

#################################################################################################################################################

# Section for HeeP physics ploting script for beam energy
elif [[ $b_flag = "true" ]]; then
    while true; do
        read -p "Do you wish to plot the ROOT files for beam energy ${BEAMENERGY}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                i=-1
                (
                echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo ""
                cd "${UTILPATH}/scripts/heep/src/"
                python3 plot_coin.py ${ROOTPREFIX_PLOT} ${BEAMENRGY} ${MAXEVENTS}
                )
                break;;
            [Nn]* )
                exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done
fi

exit 0
