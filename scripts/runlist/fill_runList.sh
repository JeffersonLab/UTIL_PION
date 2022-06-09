#!/bin/bash

# 12/08/21 - Stephen JD Kay - University of Regina
# Script to fill the run list, executed by "master" script after each run

# Set up paths depending upon location

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

# Run number and run type should be read in by the "master" script, automating the target would be more difficult, master script prompts for this
RUNNUMBER=$1
RUNTYPE=$2
TARGET=$3
if [[ ${ANATYPE} = *"Pion"* ]]; then
    RUNLIST="${UTILPATH}/runlist_pionLT_2022.csv"
elif [[ ${ANATYPE} = *"Kaon"* ]]; then
    RUNLIST="${UTILPATH}/runlist_kaonLT_2021.csv"
fi
# Need to fix paths rather than give relative paths, also need to check information is still in these files and that it can grab it correctly
KINFILE="${REPLAYPATH}/DBASE/COIN/standard.kinematics"
# Get report file based upon run type
if [[ ${RUNTYPE} = *"Prod"* ]]; then
    REPORTFILE="${REPLAYPATH}/REPORT_OUTPUT/Analysis/${ANATYPE}LT/${ANATYPE}_replay_coin_production_${RUNNUMBER}_-1.report" # Finalised
elif [[ ${RUNTYPE} = *"Lumi"* ]]; then
    REPORTFILE="${REPLAYPATH}/REPORT_OUTPUT/Analysis/Lumi/${ANATYPE}_replay_luminosity_${RUNNUMBER}_-1.report" 
elif [[ ${RUNTYPE} = *"HeePSing"* ]]; then
    REPORTFILE="${REPLAYPATH}/REPORT_OUTPUT/Analysis/HeeP/${ANATYPE}_replay_shms_production_${RUNNUMBER}_-1.report" # Finalised, all of the available info SHOULD be in the SHMS report file, don't need to look at both
elif [[ ${RUNTYPE} = *"HeePCoin"* ]]; then
    REPORTFILE="${REPLAYPATH}/REPORT_OUTPUT/Analysis/HeeP/${ANATYPE}_replay_coin_production_${RUNNUMBER}_-1.report" # Finalised 
elif [[ ${RUNTYPE} = *"fADC"* ]]; then
    REPORTFILE="${REPLAYPATH}/REPORT_OUTPUT/Analysis/${ANATYPE}LT/${ANATYPE}_replay_coin_production_${RUNNUMBER}_-1.report" # Finalised, it just uses a ${ANATYPE}LT replay
else
    REPORTFILE="${REPLAYPATH}/REPORT_OUTPUT/Analysis/General/${ANATYPE}_replay_coin_production_${RUNNUMBER}_-1.report" # CHANGE WHEN FINALISED
fi

# Get information available in standard.kinematics, execute a python script to do this for us
KINFILE_INFO=`python3 $UTILPATH/scripts/runlist/kinfile.py ${KINFILE} ${RUNNUMBER}` # The output of this python script is just a comma separated string
# Split the string we get to individual variables, easier for printing and use later
SHMS_Angle=`echo ${KINFILE_INFO} | cut -d ','  -f1` # Cut the string on , delimitter, select field (f) 1, set variable to output of command
SHMS_P=`echo ${KINFILE_INFO} | cut -d ','  -f2`
SHMS_mass=`echo ${KINFILE_INFO} | cut -d ','  -f3`
HMS_Angle=`echo ${KINFILE_INFO} | cut -d ','  -f4`
HMS_P=`echo ${KINFILE_INFO} | cut -d ','  -f5`
HMS_mass=`echo ${KINFILE_INFO} | cut -d ','  -f6`
EBeam=`echo ${KINFILE_INFO} | cut -d ','  -f7`

# Get information available in the report file
if [[ -f ${REPORTFILE} ]]; then
    if [[ ${RUNTYPE} == "Lumi" ]]; then
	REPORTFILE_INFO=`python3 $UTILPATH/scripts/runlist/reportfile_Lumi.py ${REPORTFILE}`
    elif [[ ${RUNTYPE} == "HeePSing" ]]; then
	REPORTFILE_INFO=`python3 $UTILPATH/scripts/runlist/reportfile_HeePSing.py ${REPORTFILE}`
    elif [[ ${RUNTYPE} != "HeePSing" || ${RUNTYPE} == "Lumi" ]]; then
	REPORTFILE_INFO=`python3 $UTILPATH/scripts/runlist/reportfile.py ${REPORTFILE}`
    fi
    Current=`echo ${REPORTFILE_INFO} | cut -d ',' -f1`
    PS1=`echo ${REPORTFILE_INFO} | cut -d ',' -f2`
    PS4=`echo ${REPORTFILE_INFO} | cut -d ',' -f3`
    PS5=`echo ${REPORTFILE_INFO} | cut -d ',' -f4`
    HMS_Rate=`echo ${REPORTFILE_INFO} | cut -d ',' -f5`
    SHMS_Rate=`echo ${REPORTFILE_INFO} | cut -d ',' -f6`
    COIN_Rate=`echo ${REPORTFILE_INFO} | cut -d ',' -f7`
    Charge=`echo ${REPORTFILE_INFO} | cut -d ',' -f8`
    Raw_HMS=`echo ${REPORTFILE_INFO} | cut -d ',' -f9`
    Raw_SHMS=`echo ${REPORTFILE_INFO} | cut -d ',' -f10`
    Raw_COIN=`echo ${REPORTFILE_INFO} | cut -d ',' -f11`
    EDTM=`echo ${REPORTFILE_INFO} | cut -d ',' -f12`
    Tracking=`echo ${REPORTFILE_INFO} | cut -d ',' -f13`
elif [[ ! -f ${REPORTFILE} ]]; then
    Current="ERROR"
    PS1="ERROR"
    PS4="ERROR"
    PS5="ERROR"
    HMS_Rate="ERROR"
    SHMS_Rate="ERROR"
    COIN_Rate="ERROR"
    Charge="ERROR"
    Raw_HMS="ERROR"
    Raw_SHMS="ERROR"
    Raw_COIN="ERROR"
    EDTM="ERROR"
    Tracking="ERROR"
fi

echo "========================================================================="
echo "These values autofill into the run list ..."
echo
echo "Run number: $RUNNUMBER"
echo "Run type: $RUNTYPE"
echo "Target: $TARGET"
echo "Beam energy: $EBeam"
echo "SHMS momentum: $SHMS_P"
echo "SHMS angle : $SHMS_Angle"
echo "SHMS particle mass : $SHMS_mass"
echo "HMS momentum: $HMS_P"
echo "HMS angle: $HMS_Angle"
echo "HMS particle mass : $HMS_mass"
echo "Current: $Current"
if [[ ${RUNTYPE} != "HeePSing" && ${RUNTYPE} != "Lumi" ]]; then
    echo "PS1 : $PS1"
else echo "PS2 : $PS1" # For HeepSingles we care about both ELReal triggers (2 and 4)
fi
echo "PS4 : $PS4"
echo "PS5 : $PS5"
echo "Raw HMS rate [kHz]: $HMS_Rate"
echo "Raw SHMS rate [kHz]: $SHMS_Rate"
if [[ ${RUNTYPE} != "HeePSing" && ${RUNTYPE} != "Lumi" ]]; then
    echo "Raw COIN rate [kHz]: $COIN_Rate"
else echo "Raw COIN rate [kHz] (No coin - SHMS rate duplicated): $COIN_Rate"
fi
echo "Charge [mC]: $Charge"
echo "Raw HMS: $Raw_HMS"
echo "Raw SHMS: $Raw_SHMS"
if [[ ${RUNTYPE} != "HeePSing" && ${RUNTYPE} != "Lumi" ]]; then
    echo "Raw coin: $Raw_COIN"
else echo "Raw coin (No coin - SHMS number duplicated): $Raw_COIN"
fi
echo "EDTM: $EDTM"
if [[ ${RUNTYPE} != "HeePSing" && ${RUNTYPE} != "Lumi" ]]; then
    echo "SHMS hadron tracking: $Tracking"
else echo "SHMS electron tracking: $Tracking"
fi
echo "========================================================================="

while true; do
    read -p "Do these values all look correct? (Please answer yes or no) " yn
    case $yn in
        [Yy]* ) break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
done

# Ask user for a comment
read -p "Enter number of pi/n events and/or any other comments: " Comment
Comment=$(echo "$Comment" | tr "," " ") # Remove any commas from the comment line as this will cause... issues
Comment=$(echo "$Comment" | tr ";" " ") # Remove any semicolons from the comment line as well, grammar get out!
Comment=$(echo "$Comment" | tr "\t" " ") # Tabs can go to hell too
# Need to fix widths of entries with blank space at some point, see the test file for widths (based on headers)
RUNLIST_INFO="${RUNNUMBER},${RUNTYPE},${TARGET},${EBeam},${SHMS_P},${SHMS_Angle},${HMS_P},${HMS_Angle},${Current},${PS1},${PS4},${PS5},${HMS_Rate},${SHMS_Rate},${COIN_Rate},${Charge},${Raw_HMS},${Raw_SHMS},${Raw_COIN},${EDTM},${Tracking},${Comment}"

# Check if there is already an entry for this run number, if there is, ask if you want to overwrite it, if not, print it to the file
DuplicateLines=() # Array to store line numbers of duplicated entries
LineNum=1 # Counter, starts at 1 since we skip the header
# Read run list, find any lines which already include an entry for this run number
while IFS='' read -r line || [[ -n "$line" ]]; do
    LineNum=$(($LineNum + 1 ))
    if [[ `echo ${line} | cut -d ','  -f1` == ${RUNNUMBER} ]]; then
	DuplicateLines[${#DuplicateLines[@]}]="${LineNum}"
    fi
done <  <(tail -n +2 ${RUNLIST}) # Ignores header line by using tail here

if [[ `echo "${#DuplicateLines[@]}"` != 0 ]]; then # If the array is not empty, check some stuff and check with the user if they want to remove duplicate entries
    # Ask if the user wants to remove duplicate lines, in a grammatically correct manner :)
    if [[ `echo "${#DuplicateLines[@]}"` == 1 ]]; then 
	read -p "$(echo "${#DuplicateLines[@]}") entry already found in the runlist for run ${RUNNUMBER}, delete dupliacte entry and print new entry to file? <Y/N> " prompt
    elif [[ `echo "${#DuplicateLines[@]}"` -gt 1 ]]; then
	read -p "$(echo "${#DuplicateLines[@]}") entries already found in the runlist for run ${RUNNUMBER}, delete dupliacte entries and print new entry to file? <Y/N> " prompt
    fi
    if [[ $prompt == "y" || $prompt == "Y" || $prompt == "yes" || $prompt == "Yes" ]]; then
	DeletedLines=0 # Counter to check how many lines we have deleted already
	# Loop over all line numbers identified earlier as being duplicates, delete them with a sed command
	for value in "${DuplicateLines[@]}"
	do
	    LineNum=$(($value-$DeletedLines)) # We need to account for any lines we delete as we go
	    sed -i "${LineNum}d" ${RUNLIST}
	    DeletedLines=$(($DeletedLines + 1))
	done
    echo ${RUNLIST_INFO} >> ${RUNLIST} # Print the run list info to the file
    else echo "Will not remove duplicate entries or print new entry to the file, please edit the runlist manually"
    fi
else
    echo ${RUNLIST_INFO} >> ${RUNLIST} # Print the run list info to the file
fi
