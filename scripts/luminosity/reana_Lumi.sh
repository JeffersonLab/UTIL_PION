#!/bin/bash

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
