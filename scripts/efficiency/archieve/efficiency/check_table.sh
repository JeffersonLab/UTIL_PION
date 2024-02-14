#! /bin/bash

#
# Description:
# ================================================================
# Time-stamp: "2024-01-11 14:27:53 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

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
	TIMESTMP="2024_01_11"
	ROOTPREFIX=replay_coin_heep
    else
	TIMESTMP="2024_01_11"
	ROOTPREFIX=replay_coin_production
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
    TIMESTMP="2022_07_28"
    if [[ $RUNTYPE = "HeePSing" ]]; then
	ROOTPREFIX=replay_${spec}_heep
    else
	ROOTPREFIX=replay_${spec}_production
    fi
elif [[ $s_flag = "true" && $u_flag = "true" ]]; then
    RUNTYPE=$2
    spec=$(echo "$3" | tr '[:upper:]' '[:lower:]')
    SPEC=$(echo "$spec" | tr '[:lower:]' '[:upper:]')
    RUNNUM=$4
    COLUMN=$5
    TIMESTMP="2022_07_28"
    if [[ $RUNTYPE = "HeePSing" ]]; then
	ROOTPREFIX=replay_${spec}_heep
    else
	ROOTPREFIX=replay_${spec}_production
    fi
elif [[ $u_flag = "true" ]]; then
    RUNTYPE=$2
    RUNNUM=$3
    COLUMN=$4
    if [[ $RUNTYPE = "HeePCoin" ]]; then
	TIMESTMP="2024_01_11"
	ROOTPREFIX=replay_coin_heep
    else
	TIMESTMP="2024_01_11"
	ROOTPREFIX=replay_coin_production
    fi    
else
    RUNTYPE=$1
    RUNNUM=$2
    COLUMN="All"
    if [[ $RUNTYPE = "HeePCoin" ]]; then
	TIMESTMP="2024_01_11"
	ROOTPREFIX=replay_coin_heep
    else
	TIMESTMP="2024_01_11"
	ROOTPREFIX=replay_coin_production
    fi
fi

cd "${SCRIPTPATH}/efficiency/src/"

python3 check_table.py $ROOTPREFIX $RUNTYPE $RUNNUM $TIMESTMP $COLUMN
