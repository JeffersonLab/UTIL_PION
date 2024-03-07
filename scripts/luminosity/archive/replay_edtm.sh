#!/bin/bash

echo
echo "Starting Luminosity/EDTM Replay Script"
#echo "I take as arguments the Run Number and max number of events!"
#RUNNUMBER=$1
#MAXEVENTS=$2
#MAXEVENTS=12500

#if [[ $1 -eq "" ]]; then
#    echo "I need a Run Number!"
#    exit 2
#fi

#if [[ $2 -eq "" ]]; then
#    echo "Only Run Number entered...I'll assume -1 events!" 
MAXEVENTS=-1 
#fi


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

###################################################################################################################################################
for RUNNUMBER in {14823..14824}
do
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
    
    if [ ! -f "$UTILPATH/ROOTfiles/Analysis/Lumi/${ANATYPE}_replay_luminosity_${RUNNUMBER}_${MAXEVENTS}.root" ]; then
	if [[ "${HOSTNAME}" != *"ifarm"* ]]; then
	    if [[ "${HOSTNAME}" == *"cdaq"* ]]; then
		eval "$REPLAYPATH/hcana -l -q -b \"$UTILPATH/scripts/replay/${ANATYPE}LT/replay_luminosity.C($RUNNUMBER,$MAXEVENTS)\""| tee $UTILPATH/REPORT_OUTPUT/Analysis/Lumi/${ANATYPE}_output_coin_production_Summary_${RUNNUMBER}_${MAXEVENTS}.report
	    else	
		eval "$REPLAYPATH/hcana -l -q -b \"$UTILPATH/scripts/replay/${ANATYPE}LT/replay_luminosity.C($RUNNUMBER,$MAXEVENTS)\"" 
	    fi
	elif [[ "${HOSTNAME}" == *"ifarm"* ]]; then
	    eval "$REPLAYPATH/hcana -l -q -b \"$UTILPATH/scripts/replay/${ANATYPE}LT/replay_luminosity.C($RUNNUMBER,$MAXEVENTS)\""| tee $UTILPATH/REPORT_OUTPUT/Analysis/Lumi/${ANATYPE}_output_coin_production_Summary_${RUNNUMBER}_${MAXEVENTS}.report
	fi
    else echo "Replayfile already found for this run in $UTILPATH/ROOTfiles/Analysis/Lumi/ - Skipping replay step"
    fi
done
#sleep 3

# Sets trigger windows
#echo
#echo "Running trigWindows.sh ${RUNNUMBER}..."
#cd ${UTILPATH}/scripts/trig_windows/
#source trigWindows.sh ${RUNNUMBER}
#cd ${UTILPATH}/scripts/trig_windows/
#source trigWindows.sh -p ${RUNNUMBER}

# Analyzes lumi runs
#echo
#echo "Running lumiyield.py ${RUNNUMBER} ${MAXEVENTS}..."
#cd ${UTILPATH}/scripts/luminosity/src/
#python3 lumiyield.py ${ANATYPE}_replay_luminosity ${RUNNUMBER} ${MAXEVENTS}
