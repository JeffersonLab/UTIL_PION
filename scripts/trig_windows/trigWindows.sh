#!/bin/bash

# Runs script in the ltsep python package that grabs current path enviroment
if [[ ${HOSTNAME} = *"cdaq"* ]]; then
    PATHFILE_INFO=`python3 /home/cdaq/pionLT-2021/PythonPackages3.6/lib/python3.6/site-packages/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
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

# Source stuff depending upon hostname. Change or add more as needed  
if [[ "${HOST}" = *"farm"* ]]; then
    if [[ "${HOST}" != *"ifarm"* ]]; then
	source /site/12gev_phys/softenv.sh 2.3
	source /apps/root/6.18.04/setroot_CUE.bash
    fi
    cd "$HCANAPATH"
    source "$HCANAPATH/setup.sh"
    cd "$REPLAYPATH"
    source "$REPLAYPATH/setup.sh"
elif [[ "${HOST}" = *"qcd"* ]]; then
    source "$REPLAYPATH/setup.sh" 
fi

cd "$REPLAYPATH"

# Define flags for grabbing trigger windows of one or multiple runs
while getopts 'hpa' flag; do
    case "${flag}" in
	h)
	    echo "The following flags can be called for the luminosity analysis..."
	    echo "    -h, help"
	    echo "    -p, plot trig cuts (requires additional arguments)"
	    echo "    -a, plot trig cuts for all lumi runs"
	    exit 0 ;;
	p) p_flag='true' ;;
	a) a_flag='true' ;;
	*) print_usage
	exit 1 ;;
    esac
done

# Grab inputs depending on flags used
if [[ $p_flag != "true" && $a_flag != "true" ]]; then
    echo "p I take as arguments the Run Number and max number of events!"
    RUNNUMBER=$1
    MAXEVENTS=$2
    #MAXEVENTS=12500
    if [[ $1 -eq "" ]]; then
	echo "I need a Run Number!"
	exit 2
    fi
    if [[ $2 -eq "" ]]; then
	echo "Only Run Number entered...I'll assume -1 events!" 
	MAXEVENTS=-1 
    fi
elif [[ $a_flag != "true" ]]; then
    echo "I take as arguments the Run Number and max number of events!"
    RUNNUMBER=$2
    MAXEVENTS=$3
    #MAXEVENTS=12500
    if [[ $2 -eq "" ]]; then
	echo "I need a Run Number!"
	exit 2
    fi
    if [[ $3 -eq "" ]]; then
	echo "Only Run Number entered...I'll assume -1 events!" 
	MAXEVENTS=-1 
    fi
fi
echo "Starting Luminosity Script"
# If no flags then run replays
if [[ $p_flag != "true" && $a_flag != "true" ]]; then

    ###################################################################################################################################################
    # RLT 09/24/21...Changed from 150k to full analysis. There may be issues with the currents/EDTM because the cuts may only be applying trip cuts to the first 150k events.
    # Section for luminosity replay script
    if [ ! -f "$UTILPATH/ROOTfiles/Scalers/coin_replay_scalers_${RUNNUMBER}_${MAXEVENTS}.root" ]; then
	eval "$REPLAYPATH/hcana -l -q -b \"SCRIPTS/COIN/SCALERS/replay_coin_scalers.C($RUNNUMBER,${MAXEVENTS})\""
	cd "$REPLAYPATH/CALIBRATION/bcm_current_map"
	root -b -l<<EOF 
.L ScalerCalib.C+
.x run.C("${UTILPATH}/ROOTfiles/Scalers/coin_replay_scalers_${RUNNUMBER}_${MAXEVENTS}.root")
.q  
EOF
	mv bcmcurrent_$RUNNUMBER.param $REPLAYPATH/PARAM/HMS/BCM/CALIB/bcmcurrent_$RUNNUMBER.param
	cd $REPLAYPATH
    else echo "Scaler replayfile already found for this run in $REPLAYPATH/ROOTfiles/Scalers - Skipping scaler replay step"
    fi

    sleep 3
    # SJDK 31/08/21 - Replays for luminosity analysis should output to Analysis/Lumi, for now this is probably fine
    if [ ! -f "$UTILPATH/ROOTfiles/Analysis/Lumi/${ANATYPE}_replay_luminosity_${RUNNUMBER}_${MAXEVENTS}.root" ]; then
	if [[ "${HOSTNAME}" != *"ifarm"* ]]; then
	    if [[ "${HOSTNAME}" == *"cdaq"* ]]; then
		eval "$REPLAYPATH/hcana -l -q -b \"$UTILPATH/scripts/replay/replay_luminosity.C($RUNNUMBER,$MAXEVENTS)\""| tee $UTILPATH/REPORT_OUTPUT/Analysis/Lumi/${ANATYPE}_output_coin_production_Summary_${RUNNUMBER}_${MAXEVENTS}.report
	    else	
		eval "$REPLAYPATH/hcana -l -q -b \"$UTILPATH/scripts/replay/replay_luminosity.C($RUNNUMBER,$MAXEVENTS)\"" 
	    fi
	elif [[ "${HOSTNAME}" == *"ifarm"* ]]; then 
	    eval "$REPLAYPATH/hcana -l -q -b \"$UTILPATH/scripts/replay/replay_luminosity.C($RUNNUMBER,$MAXEVENTS)\""| tee $UTILPATH/REPORT_OUTPUT/Analysis/Lumi/${ANATYPE}_output_coin_production_Summary_${RUNNUMBER}_${MAXEVENTS}.report
	fi
    else echo "Replayfile already found for this run in $UTILPATH/ROOTfiles/Analysis/Lumi/ - Skipping replay step"
    fi

    sleep 3

    cd ${UTILPATH}/scripts/trig_windows/src/
    python3 trigcuts.py Lumi ${ANATYPE}_replay_luminosity ${RUNNUMBER} ${MAXEVENTS}

fi
# Get trigger windows for a particular run
if [[ $p_flag = "true" ]]; then
    cd ${UTILPATH}/scripts/trig_windows/src/
    python3 plot_trig.py Lumi ${ANATYPE}_replay_luminosity ${RUNNUMBER} ${MAXEVENTS}
fi
# Get trigger windows for all runs
if [[ $a_flag = "true" ]]; then
    cd ${UTILPATH}/scripts/trig_windows/src/
    python3 reana_trig.py
    cd ${UTILPATH}/scripts/trig_windows/OUTPUTS/
    convert curr_${ANATYPE}_replay_* trig_${ANATYPE}_replay_* curr_trig_cuts.pdf
fi
