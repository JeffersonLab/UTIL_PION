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

while getopts 'hlus' flag; do
    case "${flag}" in
        h) 
        echo "-------------------------------------------------------------------"
        echo "./check_table.sh -{flags} {variable arguments, see help}"
        echo "-------------------------------------------------------------------"
        echo
        echo "The following flags can be called to check efficiency table..."
	echo "    If no flags called arguments are..."
	echo "        coin -> RUNTYPE=arg1 RUNNUM=arg2"
	echo "        sing -> RUNTYPE=arg1 SPEC=arg2 RUNNUM=arg3 (requires -s flag)"		
        echo "    -h, help"
        echo "    -l, list table columns"
	echo "    -u, Unique column"
	echo "        coin -> RUNTYPE=arg1 RUNNUM=arg2 COLUMN=arg3"
	echo "        sing -> RUNTYPE=arg1 SPEC=arg2 RUNNUM=arg3 COLUMN=arg4 (requires -s flag)"
	echo "    -s, single arm"
        exit 0
        ;;
	l) l_flag='true' ;;
	u) u_flag='true' ;;
	s) s_flag='true' ;;
        *) print_usage
        exit 1 ;;
    esac
done

if [[ $l_flag = "true" ]]; then
    RUNTYPE=$2
    RUNNUM=$3
    COLUMN="List"
    if [[ $RUNTYPE = "HeePCoin" ]]; then
	TIMESTMP="2024_09_17"
	ROOTPREFIX=PionLT_HeeP_coin
    elif [[ $RUNTYPE = "LumiCoin" ]]; then
        TIMESTMP="2024_09_17"
        ROOTPREFIX=PionLT_luminosity_coin
    elif [[ $RUNTYPE = "Prod" ]]; then
        TIMESTMP="2024_09_17"
        ROOTPREFIX=PionLT_coin_production
    elif [[ $RUNTYPE = "pTRIG6" ]]; then
        TIMESTMP="2024_09_17"
        ROOTPREFIX=PionLT_coin_production_pTRIG6
    else
	TIMESTMP="2024_09_17"
	ROOTPREFIX=PionLT_coin_production
    fi
    cd "${SCRIPTPATH}/efficiency/src/"
    python3 check_table.py $ROOTPREFIX $RUNTYPE $RUNNUM $TIMESTMP $COLUMN
    exit 0
fi

if [[ $s_flag = "true" ]]; then
    RUNTYPE=$2
    spec=$(echo "$3" | tr '[:upper:]' '[:lower:]')
    SPEC=$(echo "$spec" | tr '[:lower:]' '[:upper:]')
    RUNNUM=$4
    COLUMN="All"
    TIMESTMP="2024_09_17"
    if [[ $RUNTYPE = "HeePSing" ]]; then
	ROOTPREFIX=PionLT_${SPEC}_HeePSing
    elif [[ $RUNTYPE = "LumiSing" ]]; then
        ROOTPREFIX=PionLT_${SPEC}_luminosity
    else
	ROOTPREFIX=PionLT_${SPEC}_HeePSing
    fi

elif [[ $s_flag = "true" && $u_flag = "true" ]]; then
    RUNTYPE=$2
    spec=$(echo "$3" | tr '[:upper:]' '[:lower:]')
    SPEC=$(echo "$spec" | tr '[:lower:]' '[:upper:]')
    RUNNUM=$4
    COLUMN=$5
    TIMESTMP="2024_09_17"
    if [[ $RUNTYPE = "HeePSing" ]]; then
        ROOTPREFIX=PionLT_${SPEC}_HeePSing
    elif [[ $RUNTYPE = "LumiSing" ]]; then
        ROOTPREFIX=PionLT_${SPEC}_luminosity
    else
        ROOTPREFIX=PionLT_${SPEC}_HeePSing
    fi

elif [[ $u_flag = "true" ]]; then
    RUNTYPE=$2
    RUNNUM=$3
    COLUMN=$4
    if [[ $RUNTYPE = "HeePCoin" ]]; then
        TIMESTMP="2024_09_17"
        ROOTPREFIX=PionLT_HeeP_coin
    elif [[ $RUNTYPE = "LumiCoin" ]]; then
        TIMESTMP="2024_09_17"
        ROOTPREFIX=PionLT_luminosity_coin
    elif [[ $RUNTYPE = "Prod" ]]; then
        TIMESTMP="2024_09_17"
        ROOTPREFIX=PionLT_coin_production
    elif [[ $RUNTYPE = "pTRIG6" ]]; then
        TIMESTMP="2024_09_17"
        ROOTPREFIX=PionLT_coin_production_pTRIG6
    else
        TIMESTMP="2024_09_17"
        ROOTPREFIX=PionLT_coin_production
    fi    

else
    RUNTYPE=$1
    RUNNUM=$2
    COLUMN="All"
    if [[ $RUNTYPE = "HeePCoin" ]]; then
        TIMESTMP="2024_09_17"
        ROOTPREFIX=PionLT_HeeP_coin
    elif [[ $RUNTYPE = "LumiCoin" ]]; then
        TIMESTMP="2024_09_17"
        ROOTPREFIX=PionLT_luminosity_coin
    elif [[ $RUNTYPE = "Prod" ]]; then
        TIMESTMP="2024_09_17"
        ROOTPREFIX=PionLT_coin_production
    elif [[ $RUNTYPE = "pTRIG6" ]]; then
        TIMESTMP="2024_09_17"
        ROOTPREFIX=PionLT_coin_production_pTRIG6
    else
        TIMESTMP="2024_09_17"
        ROOTPREFIX=PionLT_coin_production
    fi
fi

cd "${SCRIPTPATH}/efficiency/src/"

python3 check_table.py $ROOTPREFIX $RUNTYPE $RUNNUM $TIMESTMP $COLUMN
