#! /bin/bash

#
# Description:
# ================================================================
# Time-stamp: "2024-02-09 12:10:00 junaid"
# ================================================================
#
# Author:  Muhammad Junaid <mjo147@uregina.ca>
#
# Copyright (c) junaid
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

while getopts 'hs' flag; do
    case "${flag}" in
        h) 
        echo "-------------------------------------------------------------------"
        echo "./display_table.sh -{flags} {variable arguments, see help}"
        echo "-------------------------------------------------------------------"
        echo
        echo "The following flags can be called to check efficiency table..."
	echo "    If no flags called arguments are..."
	echo "        coin -> RUNTYPE=arg1"
	echo "        sing -> RUNTYPE=arg1 SPEC=arg2 (requires -s flag)"		
        echo "    -h, help"
	echo "    -s, single arm"
        exit 0
        ;;
	s) s_flag='true' ;;
        *) print_usage
        exit 1 ;;
    esac
done

cd "${SCRIPTPATH}/efficiency/src/"

if [[ $s_flag = "true" ]]; then
    RUNTYPE=$2
    spec=$(echo "$3" | tr '[:upper:]' '[:lower:]')
    SPEC=$(echo "$spec" | tr '[:lower:]' '[:upper:]')
    TIMESTMP="2024_02_08"
    if [[ $RUNTYPE = "HeePSing" ]]; then
	ROOTPREFIX=PionLT_${SPEC}_HeePSing
        python3 plot/plot_heepsing_${SPEC}_efficiency.py ${ROOTPREFIX} ${RUNTYPE} ${TIMESTMP}
        cd "${SCRIPTPATH}/efficiency/OUTPUTS/plots"
        convert *.png "${ROOTPREFIX}_${RUNTYPE}_${TIMESTMP}.pdf"
#      evince "${RUNTYPE}_${TIMESTMP}.pdf"
        mv *.png png/
#      rm -f *.png
        exit 1
    elif [[ $RUNTYPE = "LumiSing" ]]; then
        ROOTPREFIX=PionLT_${SPEC}_luminosity
        python3 plot/plot_lumising_${SPEC}_efficiency.py ${ROOTPREFIX} ${RUNTYPE} ${TIMESTMP}
        cd "${SCRIPTPATH}/efficiency/OUTPUTS/plots"
        convert *.png "${ROOTPREFIX}_${RUNTYPE}_${TIMESTMP}.pdf"
#      evince "${RUNTYPE}_${TIMESTMP}.pdf"
        mv *.png png/
#      rm -f *.png
        exit 1
    else	
        echo "Please Provide RUNTYPE"
    fi


else
    RUNTYPE=$1
    TIMESTMP="2024_02_08"
    if [[ $RUNTYPE = "HeePCoin" ]]; then
        ROOTPREFIX=PionLT_HeeP_coin
        python3 plot/plot_heepcoin_efficiency.py ${ROOTPREFIX} ${RUNTYPE} ${TIMESTMP}
        cd "${SCRIPTPATH}/efficiency/OUTPUTS/plots"
        convert *.png "${ROOTPREFIX}_${RUNTYPE}_${TIMESTMP}.pdf"
#      evince "${RUNTYPE}_${TIMESTMP}.pdf"
        mv *.png png/
#      rm -f *.png
        exit 1
    elif [[ $RUNTYPE = "LumiCoin" ]]; then
        ROOTPREFIX=PionLT_luminosity_coin
        python3 plot/plot_lumicoin_efficiency.py ${ROOTPREFIX} ${RUNTYPE} ${TIMESTMP}
        cd "${SCRIPTPATH}/efficiency/OUTPUTS/plots"
        convert *.png "${ROOTPREFIX}_${RUNTYPE}_${TIMESTMP}.pdf"
#      evince "${RUNTYPE}_${TIMESTMP}.pdf"
        mv *.png png/
#      rm -f *.png
        exit 1
    elif [[ $RUNTYPE = "Prod" ]]; then
        ROOTPREFIX=PionLT_coin_production
        python3 plot/plot_prod_efficiency.py ${ROOTPREFIX} ${RUNTYPE} ${TIMESTMP}
        cd "${SCRIPTPATH}/efficiency/OUTPUTS/plots"
        convert *.png "${ROOTPREFIX}_${RUNTYPE}_${TIMESTMP}.pdf"
#      evince "${RUNTYPE}_${TIMESTMP}.pdf"
        mv *.png png/
#      rm -f *.png
        exit 1
    elif [[ $RUNTYPE = "pTRIG6" ]]; then
        ROOTPREFIX=PionLT_coin_production_pTRIG6
        python3 plot/plot_ptrig6_efficiency.py ${ROOTPREFIX} ${RUNTYPE} ${TIMESTMP}
        cd "${SCRIPTPATH}/efficiency/OUTPUTS/plots"
        convert *.png "${ROOTPREFIX}_${RUNTYPE}_${TIMESTMP}.pdf"
#      evince "${RUNTYPE}_${TIMESTMP}.pdf"
        mv *.png png/
#      rm -f *.png
        exit 1
    else
        echo "Please Provide RUNTYPE"
    fi
fi
