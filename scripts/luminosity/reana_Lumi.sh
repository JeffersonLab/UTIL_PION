#!/bin/bash

# Set path depending upon hostname. Change or add more as needed  
if [[ "${HOSTNAME}" = *"farm"* ]]; then  
    REPLAYPATH="/group/c-pionlt/USERS/${USER}/hallc_replay_lt"
elif [[ "${HOSTNAME}" = *"qcd"* ]]; then
    REPLAYPATH="/group/c-pionlt/USERS/${USER}/hallc_replay_lt"
elif [[ "${HOSTNAME}" = *"cdaq"* ]]; then
    REPLAYPATH="/home/cdaq/hallc-online/hallc_replay_lt"
elif [[ "${HOSTNAME}" = *"phys.uregina.ca"* ]]; then
    REPLAYPATH="/home/${USER}/work/JLab/hallc_replay_lt"
elif [[ "${HOSTNAME}" = *"trottar"* ]]; then
    REPLAYPATH="/home/trottar/Analysis/hallc_replay_lt"
fi

UTILPATH="${REPLAYPATH}/UTIL_PIONLT"

while getopts 'hrp' flag; do
    case "${flag}" in
	h)
	    echo "The following flags can be called for the luminosity analysis..."
	    echo "    -h, help"
	    echo "    -r, reanalyze all lumi runs"
	    echo "    -p, plot yield results (requires additional arguments)"
	    exit 0 ;;
	r) r_flag='true' ;;
	p) p_flag='true' ;;
	*) print_usage
	exit 1 ;;
    esac
done

PLOTINFO=$2

cd src/
if [[ $r_flag = "true" ]]; then
    echo
    echo "Reanalyzing all luminosity data..."
    LUMIFILE="${UTILPATH}/scripts/luminosity/OUTPUTS/lumi_data.csv"
    if [[ -f "${LUMIFILE}" ]]; then
	echo "Removing ${LUMIFILE}"
	rm -rf $LUMIFILE
    fi
    python3 reana_lumi.py --reana
else
    python3 reana_lumi.py
fi
if [[ $p_flag = "true" ]]; then
    if [[ ${PLOTINFO} == "" ]]; then
	echo
	echo "Need plotting arguments with -p flag!"
	exit 2
    else
	echo
	echo "Plotting yield data for ${PLOTINFO}..."
	python3 plot_yield.py ${PLOTINFO}
    fi
fi
