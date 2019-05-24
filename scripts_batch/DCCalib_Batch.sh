#!/bin/bash

### Script for running (via batch or otherwise) the DC calibration, this one script does all of the relevant steps for the calibration process
### REQUIRES two arguments, runnumber and spectrometer (HMS or SHMS, the caps are important!)
### If you want to run with LESS than all of the events, provide a third argument with # events

# The path to your hallc replay directory, change as needed
REPLAYPATH="/u/group/c-kaonlt/USERS/${USER}/hallc_replay_lt"
RUNNUMBER=$1
OPT=$2
### Check the extra folders you'll need exist, if they don't then make them

if [ ! -d "$REPLAYPATH/DBASE/COIN/HMS_DCCalib" ]; then
    mkdir "$REPLAYPATH/DBASE/COIN/HMS_DCCalib"
fi

if [ ! -d "$REPLAYPATH/DBASE/COIN/SHMS_DCCalib" ]; then
    mkdir "$REPLAYPATH/DBASE/COIN/SHMS_DCCalib"
fi

if [ ! -d "$REPLAYPATH/PARAM/HMS/DC/CALIB" ]; then
    mkdir "$REPLAYPATH/PARAM/HMS/DC/CALIB"
fi

if [ ! -d "$REPLAYPATH/PARAM/SHMS/DC/CALIB" ]; then
    mkdir "$REPLAYPATH/PARAM/SHMS/DC/CALIB"
fi

### Check you've provided the first argument
if [[ $1 -eq "" ]]; then
    echo "I need a Run Number!"
    echo "Please provide a run number as input"
    exit 2
fi

### Check you have provided the second argument correctly
if [[ ! $2 =~ ^("HMS"|"SHMS")$ ]]; then
    echo "Please specify spectrometer, HMS or SHMS"
    exit 2
fi

### Check if a third argument was provided, if not assume -1, if yes, this is max events
if [[ $3 -eq "" ]]; then
    MAXEVENTS=-1
else
    MAXEVENTS=$3
fi

if [[ $OPT == "HMS" ]]; then
    spec="hms"
    specL="h"
    elif [[ $OPT == "SHMS" ]]; then
    spec="shms"
    specL="p"
fi

# Initialize enviroment
#export OSRELEASE="Linux_CentOS7.2.1511-x86_64-gcc5.2.0"
source /site/12gev_phys/softenv.sh 2.1

# Initialize hcana, change if not running on the farm!
# Change this path if you're not on the JLab farm!
cd "/u/group/c-kaonlt/hcana/"
source "/u/group/c-kaonlt/hcana/setup.sh"
cd "$REPLAYPATH"
source "$REPLAYPATH/setup.sh"

### Run the first replay script, then, run the calibration macro
eval "$REPLAYPATH/hcana -l -q \"SCRIPTS/"$OPT"/PRODUCTION/"$OPT"DC_Calib_Coin_Pt1.C($RUNNUMBER,$MAXEVENTS)\""
sleep 30
ROOTFILE="$REPLAYPATH/ROOTfilesDCCalib/"$OPT"_DC_Calib_Pt1_"$RUNNUMBER"_"$MAXEVENTS".root" 
cd "$REPLAYPATH/CALIBRATION/dc_calib/scripts"
root -l -b -q "$REPLAYPATH/CALIBRATION/dc_calib/scripts/main_calib.C(\"$OPT\", \"$ROOTFILE\", $RUNNUMBER)"
sleep 15

### Loop checks if the new parameter files exist, returns an error if they don't
if [[ ! -f "$REPLAYPATH/CALIBRATION/dc_calib/scripts/"$OPT"_DC_cardLog_"$RUNNUMBER"/"$specL"dc_calib_"$RUNNUMBER$".param" && ! -f  "$REPLAYPATH/CALIBRATION/dc_calib/scripts/"$OPT"_DC_cardLog_"$RUNNUMBER"/"$specL"dc_tzero\
_per_wire_"$RUNNUMBER$".param" ]]; then
    echo "New parameter files not found, calibration script likely failed"
    exit 2
fi

### Copy our new parameter files to another directory
cp "$REPLAYPATH/CALIBRATION/dc_calib/scripts/"$OPT"_DC_cardLog_"$RUNNUMBER"/"$specL"dc_calib_"$RUNNUMBER$".param" "$REPLAYPATH/PARAM/"$OPT"/DC/CALIB/"$specL"dc_calib_"$RUNNUMBER$".param"
cp "$REPLAYPATH/CALIBRATION/dc_calib/scripts/"$OPT"_DC_cardLog_"$RUNNUMBER"/"$specL"dc_tzero_per_wire_"$RUNNUMBER$".param" "$REPLAYPATH/PARAM/"$OPT"/DC/CALIB/"$specL"dc_tzero_per_wire_"$RUNNUMBER$".param"
cd "$REPLAYPATH/DBASE/COIN"

### For run numbers under 6580, we use the run period 1 files as a base
### Copy these files to a new directory and rename them
### Replace info in lines 3, 37 and 38 with the path to our new files via sed commands
if [ "$RUNNUMBER" -le "6580" ]; then
    cp "$REPLAYPATH/DBASE/COIN/standard.database" "$REPLAYPATH/DBASE/COIN/"$OPT"_DCCalib/standard_$RUNNUMBER.database"
    cp "$REPLAYPATH/DBASE/COIN/general.param" "$REPLAYPATH/DBASE/COIN/"$OPT"_DCCalib/general_$RUNNUMBER.param"
    sed -i "3 s/general.param/"$OPT"_DCCalib\/general_$RUNNUMBER.param/" $REPLAYPATH"/DBASE/COIN/"$OPT"_DCCalib/standard_"$RUNNUMBER".database"
    if [[ $OPT == "HMS" ]];then
	sed -i "37 s/hdc_calib.param/CALIB\/hdc_calib_$RUNNUMBER.param/" $REPLAYPATH"/DBASE/COIN/"$OPT"_DCCalib/general_"$RUNNUMBER".param"
	sed -i "38 s/hdc_tzero_per_wire.param/CALIB\/hdc_tzero_per_wire_$RUNNUMBER.param/" $REPLAYPATH"/DBASE/COIN/"$OPT"_DCCalib/general_"$RUNNUMBER".param"
    fi
    if [[ $OPT == "SHMS" ]];then
	sed -i "71 s/dc_calib.param/CALIB\/pdc_calib_$RUNNUMBER.param/" $REPLAYPATH"/DBASE/COIN/"$OPT"_DCCalib/general_"$RUNNUMBER".param"
	sed -i "72 s/pdc_tzero_per_wire.param/CALIB\/pdc_tzero_per_wire_$RUNNUMBER.param/" $REPLAYPATH"/DBASE/COIN/"$OPT"_DCCalib/general_"$RUNNUMBER".param"
    fi
fi
### For run numbers over 6580, we use the run period 2 files as a base
### Copy these files to a new directory and rename them
### Replace info in lines 8, 37 and 38 with the path to our new files via sed commands
if [ "$RUNNUMBER" -ge "6581" ]; then
    cp "$REPLAYPATH/DBASE/COIN/standard.database" "$REPLAYPATH/DBASE/COIN/"$OPT"_DCCalib/standard_$RUNNUMBER.database"
    cp "$REPLAYPATH/DBASE/COIN/general_runperiod2.param" "$REPLAYPATH/DBASE/COIN/"$OPT"_DCCalib/general_$RUNNUMBER.param"
    sed -i "8 s/general_runperiod2.param/"$OPT"_DCCalib\/general_$RUNNUMBER.param/" $REPLAYPATH"/DBASE/COIN/"$OPT"_DCCalib/standard_"$RUNNUMBER".database"
    if [[ $OPT == "HMS" ]];then
	sed -i "37 s/hdc_calib.param/CALIB\/hdc_calib_$RUNNUMBER.param/" $REPLAYPATH"/DBASE/COIN/"$OPT"_DCCalib/general_"$RUNNUMBER".param"
	sed -i "38 s/hdc_tzero_per_wire.param/CALIB\/hdc_tzero_per_wire_$RUNNUMBER.param/" $REPLAYPATH"/DBASE/COIN/"$OPT"_DCCalib/general_"$RUNNUMBER".param"
    fi
    if [[ $OPT == "SHMS" ]]; then
	sed -i "71 s/pdc_calib.param/CALIB\/pdc_calib_$RUNNUMBER.param/" $REPLAYPATH"/DBASE/COIN/"$OPT"_DCCalib/general_"$RUNNUMBER".param"
	sed -i "72 s/pdc_tzero_per_wire.param/CALIB\/pdc_tzero_per_wire_$RUNNUMBER.param/" $REPLAYPATH"/DBASE/COIN/"$OPT"_DCCalib/general_"$RUNNUMBER".param"
    fi
fi


sleep 10
### Finally, replay again with our new parameter files
cd $REPLAYPATH
eval "$REPLAYPATH/hcana -l -q \"SCRIPTS/"$OPT"/PRODUCTION/"$OPT"DC_Calib_Coin_Pt2.C($RUNNUMBER,$MAXEVENTS)\""
