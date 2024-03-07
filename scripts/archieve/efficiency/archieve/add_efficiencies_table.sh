#! /bin/bash

#
# Description:
# ================================================================
# Time-stamp: "2022-12-30 22:25:11 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

while getopts 'hprs' flag; do
    case "${flag}" in
        h) 
        echo "-------------------------------------------------------------------"
        echo "./add_efficiencies_table.sh -{flags} {variable arguments, see help}"
        echo "-------------------------------------------------------------------"
        echo
        echo "The following flags can be called for the heep analysis..."
	echo "    If no flags called arguments are..."
	echo "        coin -> RunType=arg1"
	echo "        sing -> RunType=arg1 SPEC=arg2 (requires -s flag)"		
        echo "    -h, help"
        echo "    -p, plot efficiencies (code exits upon completion, table required)"
	echo "        coin -> RunType=arg1 DATE=arg2"
	echo "        sing -> RunType=arg1 DATE=arg2 SPEC=arg3 (requires -s flag)"
	echo "    -r, run hgcer root analysis"
	echo "        coin -> RunType=arg1"
	echo "        sing -> RunType=arg1 SPEC=arg2 (requires -s flag)"	
	echo "    -s, single arm"
        exit 0
        ;;
        p) p_flag='true' ;;
	r) r_flag='true' ;;
	s) s_flag='true' ;;
        *) print_usage
        exit 1 ;;
    esac
done

# Runs script in the ltsep python package that grabs current path enviroment
if [[ ${HOSTNAME} = *"cdaq"* ]]; then
    PATHFILE_INFO=`python3 /home/cdaq/pionLT-2021/hallc_replay_lt/UTIL_PION/bin/python/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
elif [[ "${HOSTNAME}" = *"farm"* ]]; then
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

cd "${SCRIPTPATH}/efficiency/src/"

if [[ $p_flag = "true" ]]; then
    RunType=$2
    DATE=$3
    if [[ $RunType = "HeePCoin" ]]; then
	ROOTPREFIX=replay_coin_heep
	HGCERPREFIX=${ANATYPE}_coin_replay_production
    else
	ROOTPREFIX=replay_coin_production
	HGCERPREFIX=${ANATYPE}_coin_replay_production	
    fi
    #python3 plot/plot_efficiency.py ${ROOTPREFIX} ${RunType} ${DATE}
    python3 plot/plot_efficiency_beam.py ${ROOTPREFIX} ${RunType} ${DATE}
    cd "${SCRIPTPATH}/efficiency/OUTPUTS/plots"
    convert *.png "${RunType}_${DATE}.pdf"
    evince "${RunType}_${DATE}.pdf"
    rm -f *.png
    exit 1
elif [[ $p_flag = "true" && $s_flag = "true" ]]; then
    RunType=$2
    DATE=$3
    spec=$(echo "$4" | tr '[:upper:]' '[:lower:]')
    SPEC=$(echo "$spec" | tr '[:lower:]' '[:upper:]')
    if [[ $RunType = "HeePSing" ]]; then
	ROOTPREFIX=replay_${spec}_heep
	HGCERPREFIX=${ANATYPE}_${SPEC}_replay_production
	inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/KaonLT_2018_2019/HeePSing_ALL"
	#inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/KaonLT_2018_2019/HeePSing_Test"
    else
	ROOTPREFIX=replay_${spec}_production
	HGCERPREFIX=${ANATYPE}_${SPEC}_replay_production
	inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/KaonLT_2018_2019/Prod_ALL"
	#inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/KaonLT_2018_2019/Prod_Test"
    fi
    #python3 plot/plot_efficiency.py ${ROOTPREFIX} ${RunType} ${DATE}
    python3 plot/plot_efficiency_beam.py ${ROOTPREFIX} ${RunType} ${DATE}
    cd "${SCRIPTPATH}/efficiency/OUTPUTS/plots"
    convert *.png "${RunType}_${DATE}.pdf"
    evince "${RunType}_${DATE}.pdf"
    exit 1
elif [[ $s_flag = "true" ]]; then
    RunType=$2
    spec=$(echo "$3" | tr '[:upper:]' '[:lower:]')
    SPEC=$(echo "$spec" | tr '[:lower:]' '[:upper:]')
    if [[ $RunType = "HeePSing" ]]; then
	ROOTPREFIX=replay_${spec}_heep
	HGCERPREFIX=${ANATYPE}_${SPEC}_replay_production
	inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/KaonLT_2018_2019/HeePSing_ALL"
	#inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/KaonLT_2018_2019/HeePSing_Test"
    else
	ROOTPREFIX=replay_${spec}_production
	HGCERPREFIX=${ANATYPE}_${SPEC}_replay_production	
	inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/KaonLT_2018_2019/ProductionLH2_ALL"
	#inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/KaonLT_2018_2019/Prod_Test"
    fi    
else
    RunType=$1
    if [[ $RunType = "HeePCoin" ]]; then
	ROOTPREFIX=replay_coin_heep
	HGCERPREFIX=${ANATYPE}_coin_replay_production
	inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/KaonLT_2018_2019/HeePCoin_ALL"
	#inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/KaonLT_2018_2019/HeePCoin_Test"
    else
	ROOTPREFIX=replay_coin_production
	HGCERPREFIX=${ANATYPE}_coin_replay_production	
	#inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/KaonLT_2018_2019/ProductionLH2_ALL"
	inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/KaonLT_2018_2019/Prod_ALL"
	#inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/KaonLT_2018_2019/Prod_Test"
    fi    
fi

if [[ $r_flag = "true" ]]; then
    while true; do
	read -p "Do you wish to analyse hgcer efficiencies with run list ${inputFile}? (Please answer yes or no) " yn
	case $yn in
	    [Yy]* )
		i=-1
		(
		##Reads in input file##
		while IFS='' read -r line || [[ -n "$line" ]]; do
		    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
		    echo "Run number read from file: $line"
		    echo ""
		    cd "${SCRIPTPATH}/efficiency/src/hgcer"
		    python3 hgcer.py ${HGCERPREFIX} $line -1
		done < "$inputFile"
		)
		break;;
	    [Nn]* ) 
		exit;;
	    * ) echo "Please answer yes or no.";;
	esac
    done
else
    while true; do
	read -p "Do you wish to create efficiency table for run list ${inputFile}? (Please answer yes or no) " yn
	case $yn in
	    [Yy]* )
		i=-1
		(
		##Reads in input file##
		while IFS='' read -r line || [[ -n "$line" ]]; do
		    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
		    echo "Run number read from file: $line"
		    echo ""
		    python3 efficiency_main.py ${ROOTPREFIX} $RunType $line -1
		done < "$inputFile"
		)
		break;;
	    [Nn]* ) 
		exit;;
	    * ) echo "Please answer yes or no.";;
	esac
    done
fi
