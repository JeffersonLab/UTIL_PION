#! /bin/bash

#
# Description:
# ================================================================
# Time-stamp: "2021-12-15 06:53:50 trottar"
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

#################################################################################################################################################

echo "Starting analysis of Proton events"
echo "I take as arguments the run number and number of events!"
# Input params - run number and max number of events
RUNNUMBER=$1
if [[ -z "$1" ]]; then
    echo "I need a run number"
    echo "Please provide a run number as input"
fi
MAXEVENTS=$2
if [[ -z "$2" ]]; then
    echo "Only Run Number entered...I'll assume -1 (all) events!" 
    MAXEVENTS=-1 
fi

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
	source /site/12gev_phys/softenv.sh 2.3
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

# Section for HeeP replay script
if [ ! -f "$UTILPATH/ROOTfiles/Scalers/coin_replay_scalers_${RUNNUMBER}_150000.root" ]; then
    eval "$REPLAYPATH/hcana -l -q \"SCRIPTS/COIN/SCALERS/replay_coin_scalers.C($RUNNUMBER,150000)\""
    cd "$REPLAYPATH/CALIBRATION/bcm_current_map"
    root -b -l<<EOF 
.L ScalerCalib.C+
.x run.C("${UTILPATH}/ROOTfiles/Scalers/coin_replay_scalers_${RUNNUMBER}_150000.root")
.q  
EOF
    mv bcmcurrent_$RUNNUMBER.param $REPLAYPATH/PARAM/HMS/BCM/CALIB/bcmcurrent_$RUNNUMBER.param
    cd $REPLAYPATH

    echo ""                                                                                                                                                                           
    echo ""                                                                                                                                                                                   
    echo ""                                                                                                                                            
    echo "-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|"                           
    echo ""                                                   
    echo "     _      "
    echo "   _| |     "
    echo " _| | |     "
    echo "| | | |     "
    echo "| | | | __  "
    echo "| | | |/  \ "
    echo "|       /\ \\"
    echo "|      /  \/ Good soup!"
    echo "|      \  /\\"
    echo "|       \/ /"
    echo " \        / "
    echo "  |     /   "
    echo "  |    |    "                          
    echo ""          
    echo "-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|"                                               
    echo ""                                                                                                                                                
    echo ""                                                                                    
    echo ""                         

else echo "Scaler replayfile already found for this run in $REPLAYPATH/ROOTfiles/Scalers - Skipping scaler replay step"
fi

sleep 3

# SJDK 31/08/21 - Replays for HeeP analysis should output to Analysis/HeeP, for now this is probably fine
if [ ! -f "$UTILPATH/ROOTfiles/Analysis/HeeP/${ANATYPE}_coin_replay_production_${RUNNUMBER}_${MAXEVENTS}.root" ]; then
    if [[ "${HOSTNAME}" != *"ifarm"* ]]; then
	if [[ "${HOSTNAME}" == *"cdaq"* ]]; then
	    eval "$REPLAYPATH/hcana -l -q \"$UTILPATH/scripts/replay/${ANATYPE}LT/replay_coin_heep.C($RUNNUMBER,$MAXEVENTS)\""| tee $UTILPATH/REPORT_OUTPUT/Analysis/HeeP/${ANATYPE}_output_coin_production_Summary_${RUNNUMBER}_${MAXEVENTS}.report
	else	
	    eval "$REPLAYPATH/hcana -l -q \"$UTILPATH/scripts/replay/${ANATYPE}LT/replay_coin_heep.C($RUNNUMBER,$MAXEVENTS)\"" 
	fi
    elif [[ "${HOSTNAME}" == *"ifarm"* ]]; then
	eval "$REPLAYPATH/hcana -l -q \"$UTILPATH/scripts/replay/${ANATYPE}LT/replay_coin_heep.C($RUNNUMBER,$MAXEVENTS)\""| tee $UTILPATH/REPORT_OUTPUT/Analysis/HeeP/${ANATYPE}_output_coin_production_Summary_${RUNNUMBER}_${MAXEVENTS}.report
    fi
else echo "Replayfile already found for this run in $UTILPATH/ROOTfiles/Analysis/HeeP/ - Skipping replay step"
fi

sleep 3

################################################################################################################################                                                                                   
# Section for HeeP analysis script
if [ -f "${UTILPATH}/OUTPUT/Analysis/HeeP/${RUNNUMBER}_${MAXEVENTS}_Analysed_Data.root" ]; then
    read -p "${ANATYPE} production analysed file already exists, do you want to reprocess it? <Y/N> " option1
    if [[ $option1 == "y" || $option1 == "Y" || $option1 == "yes" || $option1 == "Yes" ]]; then
	rm "${UTILPATH}/OUTPUT/Analysis/HeeP/${RUNNUMBER}_${MAXEVENTS}_Analysed_Data.root"
	echo "Reprocessing"
	python3 ${UTILPATH}/scripts/heep/src/coinyield.py ${ANATYPE}_coin_replay_production ${RUNNUMBER} ${MAXEVENTS}
    else
	echo "Skipping python analysis script step"
    fi
elif [ ! -f "${UTILPATH}/OUTPUT/Analysis/HeeP/${RUNNUMBER}_${MAXEVENTS}_Analysed_Data.root" ]; then
	python3 ${UTILPATH}/scripts/heep/src/coinyield.py ${ANATYPE}_coin_replay_production ${RUNNUMBER} ${MAXEVENTS}
else echo "Analysed root file already found in ${UTILPATH}/OUTPUT/Analysis/HeeP/ - Skipped python analyzer script step"
fi

sleep 3

##################################################################################################################################

# Section for HeeP physics ploting script
# 23/09/21 - SJDK - Changed the ordering of the arguments given to the python script to make them consistent
if [ -f "${UTILPATH}/OUTPUT/Analysis/HeeP/${RUNNUMBER}_${MAXEVENTS}_Output_Data.root" ]; then
    read -p "HeeP output plots file already exists, do you want to reprocess it? <Y/N> " option2
    if [[ $option2 == "y" || $option2 == "Y" || $option2 == "yes" || $option2 == "Yes" ]]; then
	rm "${UTILPATH}/OUTPUT/Analysis/HeeP/${RUNNUMBER}_${MAXEVENTS}_Output_Data.root"
	echo "Reprocessing"
	python3 ${UTILPATH}/scripts/heep/src/plot_coin.py Analysed_Data ${RUNNUMBER} ${MAXEVENTS} 
    else
	echo "Skipping python HeeP plotting script step"
    fi
elif [ ! -f  "${UTILPATH}/OUTPUT/Analysis/HeeP/${RUNNUMBER}_${MAXEVENTS}_Output_Data.root" ]; then
	python3 ${UTILPATH}/scripts/heep/src/plot_coin.py Analysed_Data ${RUNNUMBER} ${MAXEVENTS} 
else echo "${ANATYPE} physics output root file already found in ${UTILPATH}/OUTPUT/Analysis/HeeP/ - Skipped python plotting script step"
fi
evince "${UTILPATH}/OUTPUT/Analysis/HeeP/${RUNNUMBER}_${MAXEVENTS}_sw_heep_Proton_Analysis_Distributions.pdf" &
exit 0
