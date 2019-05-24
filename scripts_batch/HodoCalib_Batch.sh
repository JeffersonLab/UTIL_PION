#!/bin/bash

### Script for running (via batch or otherwise) thehodoscope calibration, this one script does all of the relevant steps for the calibration proces
### Note that the second part also has an additional bit where it checks for a database file based upon the run number

# The path to your hallc replay directory, change as needed                                                                                                                                                        
REPLAYPATH="/u/group/c-kaonlt/USERS/${USER}/hallc_replay_lt"
RUNNUMBER=$1
OPT=$2
### Check the extra folders you'll need exist, if they don't then make them                                                                                                                                       
if [ ! -d "$REPLAYPATH/DBASE/COIN/HMS_HodoCalib" ]; then
    mkdir "$REPLAYPATH/DBASE/COIN/HMS_HodoCalib"
fi

if [ ! -d "$REPLAYPATH/DBASE/COIN/SHMS_HodoCalib" ]; then
    mkdir "$REPLAYPATH/DBASE/COIN/SHMS_HodoCalib"
fi

if [ ! -d "$REPLAYPATH/PARAM/HMS/HODO/Calibration" ]; then
    mkdir "$REPLAYPATH/PARAM/HMS/HODO/Calibration"
fi

if [ ! -d "$REPLAYPATH/PARAM/SHMS/HODO/Calibration" ]; then
    mkdir "$REPLAYPATH/PARAM/SHMS/HODO/Calibration"
fi

if [ ! -d "$REPLAYPATH/CALIBRATION/hms_hodo_calib/Calibration_Plots" ]; then
    mkdir "$REPLAYPATH/CALIBRATION/hms_hodo_calib/Calibration_Plots"
fi

if [ ! -d "$REPLAYPATH/CALIBRATION/shms_hodo_calib/Calibration_Plots" ]; then
    mkdir "$REPLAYPATH/CALIBRATION/shms_hodo_calib/Calibration_Plots"
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

#Initialize enviroment
#export OSRELEASE="Linux_CentOS7.2.1511-x86_64-gcc5.2.0"
source /site/12gev_phys/softenv.sh 2.1

# Initialize hcana, change the path if not running on the farm!
cd "/u/group/c-kaonlt/hcana/"
source "/u/group/c-kaonlt/hcana/setup.sh"
cd "$REPLAYPATH"
source "$REPLAYPATH/setup.sh"

eval "$REPLAYPATH/hcana -l -q \"SCRIPTS/"$OPT"/PRODUCTION/"$OPT"Hodo_Calib_Coin_Pt1.C($RUNNUMBER,$MAXEVENTS)\""
sleep 30
ROOTFILE="$REPLAYPATH/ROOTfilesHodoCalib/"$OPT"_Hodo_Calib_Pt1_"$RUNNUMBER"_"$MAXEVENTS".root" 

if [[ $OPT == "HMS" ]]; then
    spec="hms"
    specL="h"
    elif [[ $OPT == "SHMS" ]]; then
    spec="shms"
    specL="p"
fi

cd "$REPLAYPATH/CALIBRATION/"$spec"_hodo_calib/"
sleep 30
root -l -q -b "$REPLAYPATH/CALIBRATION/"$spec"_hodo_calib/timeWalkHistos.C(\"$ROOTFILE\", $RUNNUMBER, \"$spec\")"
sleep 60
root -l -q -b "$REPLAYPATH/CALIBRATION/"$spec"_hodo_calib/timeWalkCalib.C($RUNNUMBER)"
sleep 30

# After executing first two root scripts, should have a new .param file so long as scripts ran ok, IF NOT THEN EXIT
if [ ! -f "$REPLAYPATH/PARAM/"$OPT"/HODO/"$specL"hodo_TWcalib_$RUNNUMBER.param" ]; then
    echo ""$specL"hodo_TWCalib_$RUNNUMBER.param not found, calibration script likely failed"
    exit 2
fi

### Now we set up the second replay by making some .database and .param files for them
cd "$REPLAYPATH/DBASE/COIN"
if [ "$RUNNUMBER" -le "6580" ]; then
    # Copy our normal ones
    cp "$REPLAYPATH/DBASE/COIN/standard.database" "$REPLAYPATH/DBASE/COIN/"$OPT"_HodoCalib/standard_$RUNNUMBER.database"
    cp "$REPLAYPATH/DBASE/COIN/general.param" "$REPLAYPATH/DBASE/COIN/"$OPT"_HodoCalib/general_$RUNNUMBER.param"
    # Use sed to replace the strings, 3 means line 3, note for sed to work with a variable we need to use "", \ is an ignore character which we need to get the \ in there syntax "Line# s/TEXT TO REPLACE/REPLACEMENT/" FILE
    sed -i "3 s/general.param/"$OPT"_HodoCalib\/general_$RUNNUMBER.param/" $REPLAYPATH"/DBASE/COIN/"$OPT"_HodoCalib/standard_$RUNNUMBER.database"
    if [[ $OPT == "HMS" ]]; then
	sed -i "40 s/hhodo_TWcalib.param/hhodo_TWcalib_$RUNNUMBER.param/" $REPLAYPATH/DBASE/COIN/HMS_HodoCalib/general_$RUNNUMBER.param 
    elif [[ $OPT == "SHMS" ]]; then
	sed -i "74 s/phodo_TWcalib.param/phodo_TWcalib_$RUNNUMBER.param/" $REPLAYPATH/DBASE/COIN/SHMS_HodoCalib/general_$RUNNUMBER.param
    fi
fi

if [ "$RUNNUMBER" -ge "6581" ]; then
    cp "$REPLAYPATH/DBASE/COIN/standard.database" "$REPLAYPATH/DBASE/COIN/"$OPT"_HodoCalib/standard_$RUNNUMBER.database"
    cp "$REPLAYPATH/DBASE/COIN/general_runperiod2.param" "$REPLAYPATH/DBASE/COIN/"$OPT"_HodoCalib/general_$RUNNUMBER.param"
    sed -i "8 s/general_runperiod2.param/"$OPT"_HodoCalib\/general_$RUNNUMBER.param/" $REPLAYPATH"/DBASE/COIN/"$OPT"_HodoCalib/standard_$RUNNUMBER.database"
    if [[ $OPT == "HMS" ]]; then
	sed -i "40 s/hhodo_TWcalib.param/hhodo_TWcalib_$RUNNUMBER.param/" $REPLAYPATH/DBASE/COIN/HMS_HodoCalib/general_$RUNNUMBER.param 
    elif [[ $OPT == "SHMS" ]]; then
	sed -i "74 s/phodo_TWcalib.param/phodo_TWcalib_$RUNNUMBER.param/" $REPLAYPATH/DBASE/COIN/SHMS_HodoCalib/general_$RUNNUMBER.param
    fi
fi

# Back to the main directory
cd "$REPLAYPATH"                                
# Off we go again replaying
sleep 30
eval "$REPLAYPATH/hcana -l -q \"SCRIPTS/"$OPT"/PRODUCTION/"$OPT"Hodo_Calib_Coin_Pt2.C($RUNNUMBER,$MAXEVENTS)\""
sleep 30
# Clean up the directories of our generated files
mv "$REPLAYPATH/PARAM/"$OPT"/HODO/"$specL"hodo_TWcalib_$RUNNUMBER.param" "$REPLAYPATH/PARAM/"$OPT"/HODO/Calibration/"$specL"hodo_TWcalib_$RUNNUMBER.param"
mv "$REPLAYPATH/CALIBRATION/"$spec"_hodo_calib/timeWalkHistos_"$RUNNUMBER".root" "$REPLAYPATH/CALIBRATION/"$spec"_hodo_calib/Calibration_Plots/timeWalkHistos_"$RUNNUMBER".root"

cd "$REPLAYPATH/CALIBRATION/"$spec"_hodo_calib/"
# Define the path to the second replay root file
ROOTFILE2="$REPLAYPATH/ROOTfilesHodoCalib/"$OPT"_Hodo_Calib_Pt2_"$RUNNUMBER"_"$MAXEVENTS".root"
# Execute final script
sleep 30
root -l -q -b "$REPLAYPATH/CALIBRATION/"$spec"_hodo_calib/fitHodoCalib.C(\"$ROOTFILE2\", $RUNNUMBER)" 
sleep 30
# Check our new file exists, if not exit, if yes, move it
if [ ! -f "$REPLAYPATH/PARAM/"$OPT"/HODO/"$specL"hodo_Vpcalib_$RUNNUMBER.param" ]; then
    echo ""$specL"hodo_Vpcalib_$RUNNUMBER.param not found, calibration script likely failed"
    exit 2
fi

mv "$REPLAYPATH/PARAM/"$OPT"/HODO/"$specL"hodo_Vpcalib_$RUNNUMBER.param" "$REPLAYPATH/PARAM/"$OPT"/HODO/Calibration/"$specL"hodo_Vpcalib_$RUNNUMBER.param"

# Check our new file exists, if not exit, if yes, move it
if [ ! -f "$REPLAYPATH/CALIBRATION/"$spec"_hodo_calib/HodoCalibPlots_$RUNNUMBER.root" ]; then
    echo "HodoCalibPlots_$RUNNUMBER.root not found, calibration script likely failed"
    exit 2
fi

mv "$REPLAYPATH/CALIBRATION/"$spec"_hodo_calib/HodoCalibPlots_$RUNNUMBER.root" "$REPLAYPATH/CALIBRATION/"$spec"_hodo_calib/Calibration_Plots/HodoCalibPlots_$RUNNUMBER.root"

### Now we set up the third replay by edditing our general.param file
cd "$REPLAYPATH/DBASE/COIN"

if [[ $OPT == "HMS" ]]; then
    sed -i "40 s/hhodo_TWcalib_$RUNNUMBER.param/Calibration\/hhodo_TWcalib_$RUNNUMBER.param/" $REPLAYPATH/DBASE/COIN/HMS_HodoCalib/general_$RUNNUMBER.param 
    sed -i "41 s/hhodo_Vpcalib.param/Calibration\/hhodo_Vpcalib_$RUNNUMBER.param/" $REPLAYPATH/DBASE/COIN/HMS_HodoCalib/general_$RUNNUMBER.param
elif [[ $OPT == "SHMS" ]]; then	
    sed -i "74 s/phodo_TWcalib_$RUNNUMBER.param/Calibration\/phodo_TWcalib_$RUNNUMBER.param/" $REPLAYPATH/DBASE/COIN/SHMS_HodoCalib/general_$RUNNUMBER.param 
    sed -i "75 s/phodo_Vpcalib.param/Calibration\/phodo_Vpcalib_$RUNNUMBER.param/" $REPLAYPATH/DBASE/COIN/SHMS_HodoCalib/general_$RUNNUMBER.param
fi

sleep 30
cd "$REPLAYPATH"
eval "$REPLAYPATH/hcana -l -q \"SCRIPTS/"$OPT"/PRODUCTION/"$OPT"Hodo_Calib_Coin_Pt3.C($RUNNUMBER,$MAXEVENTS)\""
sleep 30
