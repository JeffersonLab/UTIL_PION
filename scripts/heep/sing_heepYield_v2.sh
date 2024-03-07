#! /bin/bash

#
# Description:
# ================================================================
# Time-stamp: "2021-12-15 06:53:57 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
##################################################################################
# Created - 10/July/2021, Author - Muhammad Junaid, University of Regina, Canada
##################################################################################
# This version of script is for shift workers at JLAB
# Executes the replay script and python analysis script and at the end python plotting script
# To run this script, execute ./scriptname $RUNNUMBER$
# Version 2 - 20/05/22 - Edited out scaler replay and BCM "calibration" SJDK. Removed soup.
#################################################################################################################################################

echo "Starting analysis of Proton events"
echo "I take as arguments the run number and number of events!"
# Input params - run number and max number of events
spec=$1
if [[ -z "$1" ]]; then
    echo "Please specify which spectrometer (hms or shms)" 
fi
RUNNUMBER=$2
if [[ -z "$2" ]]; then
    echo "I need a run number"
    echo "Please provide a run number as input"
fi
MAXEVENTS=$3
if [[ -z "$3" ]]; then
    echo "Only Run Number entered...I'll assume -1 (all) events!" 
    MAXEVENTS=-1 
fi

SPEC=$(echo "$spec" | tr '[:lower:]' '[:upper:]')

# Runs script in the ltsep python package that grabs current path enviroment
if [[ ${HOSTNAME} = *"cdaq"* ]]; then
    PATHFILE_INFO=`python3 /home/cdaq/pionLT-2021/hallc_replay_lt/UTIL_PION/bin/python/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
elif [[ "${HOSTNAME}" = *"farm"* ]]; then
    PATHFILE_INFO=`python3 /u/home/${USER}/.local/lib/python3.4/site-packages/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
elif [[ "${HOSTNAME}" = *"qcd"* ]]; then
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

#################################################################################################################################################

# Set path depending upon hostname. Change or add more as needed  
if [[ "${HOSTNAME}" = *"farm"* ]]; then  
    if [[ "${HOSTNAME}" != *"ifarm"* ]]; then
	#source /site/12gev_phys/softenv.sh 2.3
	source /apps/root/6.18.04/setroot_CUE.bash
    fi
    cd "$HCANAPATH"
    source "$HCANAPATH/setup.sh"
    cd "$REPLAYPATH"
    source "$REPLAYPATH/setup.sh"
elif [[ "${HOSTNAME}" = *"qcd"* ]]; then
    cd "$REPLAYPATH"
    source "$REPLAYPATH/setup.sh" 
fi

cd $REPLAYPATH

###################################################################################################################################################
# Section for pion replay script

if [ ! -f "$UTILPATH/ROOTfiles/Analysis/HeeP/${ANATYPE}_${spec}_replay_production_${RUNNUMBER}_${MAXEVENTS}.root" ]; then
    if [[ "${HOSTNAME}" != *"ifarm"* ]]; then
	if [[ "${HOSTNAME}" == *"cdaq"* ]]; then
	    eval "$REPLAYPATH/hcana -l -q \"$UTILPATH/scripts/replay/${ANATYPE}LT/replay_${spec}_heep.C($RUNNUMBER,$MAXEVENTS)\""| tee $UTILPATH/REPORT_OUTPUT/Analysis/HeeP/${ANATYPE}_output_${spec}_production_Summary_${RUNNUMBER}_${MAXEVENTS}.report
	else	
	    eval "$REPLAYPATH/hcana -l -q \"$UTILPATH/scripts/replay/${ANATYPE}LT/replay_${spec}_heep.C($RUNNUMBER,$MAXEVENTS)\"" 
	fi
    elif [[ "${HOSTNAME}" == *"ifarm"* ]]; then
	eval "$REPLAYPATH/hcana -l -q \"$UTILPATH/scripts/replay/${ANATYPE}LT/replay_${spec}_heep.C($RUNNUMBER,$MAXEVENTS)\""| tee $UTILPATH/REPORT_OUTPUT/Analysis/HeeP/${ANATYPE}_output_${spec}_production_Summary_${RUNNUMBER}_${MAXEVENTS}.report
    fi
else echo "Replayfile already found for this run in $UTILPATH/ROOTfiles/Analysis/HeeP/ - Skipping replay step"
fi

sleep 3

################################################################################################################################                                                                                   
# Section for pion analysis script
if [ -f "${UTILPATH}/OUTPUT/Analysis/HeeP/${SPEC}_${RUNNUMBER}_${MAXEVENTS}_Analysed_Data.root" ]; then
    read -p "${SPEC} HeeP Singles production analyzed file already exits, do you want to reprocess it? <Y/N> " option1
    if [[ $option1 == "y" || $option1 == "Y" || $option1 == "yes" || $option1 == "Yes" ]]; then
	rm "${UTILPATH}/OUTPUT/Analysis/HeeP/${SPEC}_${RUNNUMBER}_${MAXEVENTS}_Analysed_Data.root"
	echo "Reprocessing"
	python3 ${UTILPATH}/scripts/heep/src/singyield.py ${ANATYPE}_${spec}_replay_production ${RUNNUMBER} ${MAXEVENTS} ${SPEC}
    else
	echo "Skipping python analysis script step"
    fi
elif [ ! -f "${UTILPATH}/OUTPUT/Analysis/HeeP/${SPEC}_${RUNNUMBER}_${MAXEVENTS}_Analysed_Data.root" ]; then
	python3 ${UTILPATH}/scripts/heep/src/singyield.py ${ANATYPE}_${spec}_replay_production ${RUNNUMBER} ${MAXEVENTS} ${SPEC}
else echo "Analysed root file already found in ${UTILPATH}/OUTPUT/Analysis/HeeP/ - Skipped python analyzer script step"
fi

sleep 3

##################################################################################################################################
# Section for pion physics ploting script
# 23/09/21 - SJDK - Changed the ordering of the arguments given to the python script to make them consistent
if [ -f "${UTILPATH}/OUTPUT/Analysis/HeeP/${SPEC}_${RUNNUMBER}_${MAXEVENTS}_Output_Data.root" ]; then
    read -p "${SPEC} HeeP Singles physics output file already exits, you want to reprocess it? <Y/N> " option2
    if [[ $option2 == "y" || $option2 == "Y" || $option2 == "yes" || $option2 == "Yes" ]]; then
	rm "${UTILPATH}/OUTPUT/Analysis/HeeP/${SPEC}_${RUNNUMBER}_${MAXEVENTS}_Output_Data.root"
	echo "Reprocessing"
	python3 ${UTILPATH}/scripts/heep/src/plot_sing.py Analysed_Data ${RUNNUMBER} ${MAXEVENTS} ${SPEC}
    else
	echo "Skipping python physics plotting script step"
    fi
elif [ ! -f  "${UTILPATH}/OUTPUT/Analysis/HeeP/${SPEC}_${RUNNUMBER}_${MAXEVENTS}_Output_Data.root" ]; then
	python3 ${UTILPATH}/scripts/heep/src/plot_sing.py Analysed_Data ${RUNNUMBER} ${MAXEVENTS} ${SPEC}
else echo "${ANATYPE} physics output root file already found in ${UTILPATH}/OUTPUT/Analysis/HeeP/ - Skipped python output script step"
fi
evince "${UTILPATH}/OUTPUT/Analysis/HeeP/${RUNNUMBER}_${MAXEVENTS}_sw_heep_${SPEC}_Analysis_Distributions.pdf" &
exit 0
