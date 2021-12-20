#!/bin/bash

echo "Starting Particle ID Script"
echo "I take as arguments the Run Number and max number of events!"
RUNNUMBER=$1
MAXEVENTS=$2
# MAXEVENTS=50000

if [[ $1 -eq "" ]]; then
    echo "I need a Run Number!"
    exit 2
fi

if [[ $2 -eq "" ]]; then
    echo "Only Run Number entered...I'll assume -1 events!" 
    MAXEVENTS=-1 
fi


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


# Set path depending upon hostname. Change or add more as needed  
if [[ "${HOSTNAME}" = *"farm"* ]]; then  
    if [[ "${HOSTNAME}" != *"ifarm"* ]]; then
	source /site/12gev_phys/softenv.sh 2.3
    fi
    cd "$HCANAPATH"
    source "$HCANAPATH/hcana/setup.sh"
    cd "$REPLAYPATH"
    source "$REPLAYPATH/setup.sh"
elif [[ "${HOSTNAME}" = *"qcd"* ]]; then
    source /site/12gev_phys/softenv.sh 2.3
    cd "$HCANAPATH"
    source "$HCANAPATH/hcana/setup.sh"
    cd "$REPLAYPATH"
    source "$REPLAYPATH/setup.sh" 
fi

cd ${REPLAYPATH}

echo -e "\n\nStarting Replay Script\n\n"
./hcana -q "${UTILPATH}/scripts/pid/src/replay/replay_pid_coin_offline.C($RUNNUMBER,$MAXEVENTS)"

source /apps/root/6.18.04/setroot_CUE.bash
cd ${UTILPATH}/scripts/pid/src/
python3 pid_eff.py ${RUNNUMBER} ${MAXEVENTS}

#cd ${UTILPATH}/UTIL_PION/scripts/pid/OUTPUTS/
#convert noID_hms_cer_${RUNNUMBER}.png PID_hms_cer_${RUNNUMBER}.png noID_hms_cal_${RUNNUMBER}.png PID_hms_cal_${RUNNUMBER}.png noID_shms_hgcer_${RUNNUMBER}.png PID_shms_hgcer_${RUNNUMBER}.png noID_shms_aero_${RUNNUMBER}.png PID_shms_aero_${RUNNUMBER}.png noID_shms_cal_${RUNNUMBER}.png PID_shms_cal_${RUNNUMBER}.png pid_plots_${RUNNUMBER}.pdf
#rm -rf *.png

cd ${UTILPATH}/scripts/pid/src/
python3 csv2root.py

