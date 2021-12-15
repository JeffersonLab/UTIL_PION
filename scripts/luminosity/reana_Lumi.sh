#!/bin/bash

# Runs script in the ltsep python package that grabs current path enviroment
if [[ ${HOSTNAME} = *"cdaq"* ]]; then
    PATHFILE_INFO=`python3 /home/cdaq/pionLT-2021/hallc_replay_lt/UTIL_PION/bin/python/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
elif [[ "${HOSTNAME}" = *"farm"* ]]; then
    PATHFILE_INFO=`python3 /u/home/${USER}/.local/lib/python3.4/site-packages/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
fi

# Split the string we get to individual variables, easier for printing and use later
HCANAPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f1` # Cut the string on , delimitter, select field (f) 1, set variable to output of command
REPLAYPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f2`
UTILPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f3`
PACKAGEPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f4`
OUTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f5`
ROOTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f6`
REPORTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f7`
CUTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f8`
PARAMPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f9`
SCRIPTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f10`
ANATYPE=`echo ${PATHFILE_INFO} | cut -d ','  -f11`
USER=`echo ${PATHFILE_INFO} | cut -d ','  -f12`
HOST=`echo ${PATHFILE_INFO} | cut -d ','  -f13`

# Flags for plotting yield or reanalyzing all data
while getopts 'hrpy' flag; do
    case "${flag}" in
	h)
	    echo "The following flags can be called for the luminosity analysis..."
	    echo "    -h, help"
	    echo "    -r, reanalyze all lumi runs"
	    echo "    -p, plot yield results (requires additional arguments)"
	    echo "    -y, plot pid to see cuts"
	    exit 0 ;;
	r) r_flag='true' ;;
	p) p_flag='true' ;;
	y) y_flag='true' ;;
	*) print_usage
	exit 1 ;;
    esac
done

PLOTINFO=$2

cd src/
if [[ $r_flag = "true" ]]; then
    echo
    echo "Reanalyzing all luminosity data..."
    python3 reana_lumi.py --reana
else
    python3 reana_lumi.py
fi
if [[ $y_flag = "true" ]]; then
    if [[ ${PLOTINFO} == "" ]]; then
	echo
	echo "Need plotting arguments with -y flag!"
	exit 2
    else
	echo
	echo "Plotting yield data for ${PLOTINFO}..."
	python3 plot_yield.py ${PLOTINFO}
    fi
fi
if [[ $p_flag = "true" ]]; then
    source /site/12gev_phys/softenv.sh 2.3
    source /apps/root/6.18.04/setroot_CUE.bash
    if [[ ${PLOTINFO} == "" ]]; then
	echo
	echo "Need run number with -p flag!"
	exit 2
    else
	echo
	echo "Plotting yield data for ${PLOTINFO}..."
	python3 plot_pid.py ${ANATYPE}_replay_luminosity ${PLOTINFO} -1
    fi
fi
