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

# Runs script in the ltsep python package that grabs current path enviroment
if [[ ${HOSTNAME} = *"cdaq"* ]]; then
    PATHFILE_INFO=`python3 /home/cdaq/pionLT-2021/hallc_replay_lt/UTIL_PION/bin/python/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
elif [[ "${HOSTNAME}" = *"farm"* ]]; then
#    PATHFILE_INFO=`python3 /u/home/${USER}/.local/lib/python3.4/site-packages/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
    PATHFILE_INFO=`python3 /u/group/c-pionlt/USERS/${USER}/replay_lt_env/lib/python3.9/site-packages/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string

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

while getopts 'hsc' flag; do
    case "${flag}" in
        h) 
        echo "-------------------------------------------------------------------"
        echo "./display_table.sh -{flags} {variable arguments, see help}"
        echo "-------------------------------------------------------------------"
        echo
        echo "The following flags can be called to check efficiency table..."
	echo "    If no flags called arguments are..."
	echo "    -c  coin -> RUNTYPE=arg1 RunList = arg2 (requires -c flag)"
        echo "    -h, help"
	echo "    -s, sing -> RUNTYPE=arg1 RunList = arg2 SPEC=arg3 (requires -s flag)"
        exit 0
        ;;
	s) s_flag='true' ;;
        c) c_flag='true' ;;
        *) print_usage
        exit 1 ;;
    esac
done

RunList=$3

cd "${SCRIPTPATH}/efficiency/src/"

if [[ $s_flag = "true" ]]; then
    RUNTYPE=$2
    spec=$(echo "$4" | tr '[:upper:]' '[:lower:]')
    SPEC=$(echo "$spec" | tr '[:lower:]' '[:upper:]')
    TIMESTMP="2024_05_28"
    if [[ $RUNTYPE = "HeePSing" ]]; then
	ROOTPREFIX=PionLT_${SPEC}_HeePSing
        inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
    elif [[ $RUNTYPE = "LumiSing" ]]; then
        ROOTPREFIX=PionLT_${SPEC}_Lumi
        inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
    else	
        echo "Please Provide RUNTYPE"
    fi

elif [[ $c_flag = "true" ]]; then
    RUNTYPE=$2
    TIMESTMP="2024_05_28"
    if [[ $RUNTYPE = "HeePCoin" ]]; then
        ROOTPREFIX=PionLT_HeeP_coin
        inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
    elif [[ $RUNTYPE = "LumiCoin" ]]; then
        ROOTPREFIX=PionLT_Lumi_coin
        inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
    elif [[ $RUNTYPE = "Prod" ]]; then
        ROOTPREFIX=PionLT_coin_production
        inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
    elif [[ $RUNTYPE = "pTRIG6" ]]; then
        ROOTPREFIX=PionLT_coin_production_pTRIG6
        inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
    else
        echo "Please Provide RUNTYPE"
    fi
else
    echo "Please Provide RUNTYPE, RunList."
fi

while true; do
      read -p "Do you wish to create efficiency plot for run list ${inputFile}? (Please answer yes or no) " yn
      case $yn in
           [Yy]* )
               i=-1
               (
               ##Reads in input file##
               while IFS='' read -r line || [[ -n "$line" ]]; do
                   echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                   echo "Run number read from file: $line"
                   echo ""
                   python3 plot/plot_heepcoin_detector_efficiency.py ${ROOTPREFIX} ${RUNTYPE} ${TIMESTMP} $line -1
               done < "$inputFile"
               )
               break;;
           [Nn]* )
               exit;;
           * ) echo "Please answer yes or no.";;
      esac
done

#        convert *.png "${ROOTPREFIX}_${RUNTYPE}_${TIMESTMP}.pdf"
#      evince "${RUNTYPE}_${TIMESTMP}.pdf"
#        mv *.png png/
#      rm -f *.png
exit 1
