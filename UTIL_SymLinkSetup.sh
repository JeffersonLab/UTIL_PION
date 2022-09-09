#!/bin/bash

# 10/08/21 - Stephen JD Kay - University of Regina
# This script sets up the relevant sym links for running analysis scripts
# This version is specific for cdaq


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

# Make the sym links
echo "Sym links will now be made for - "
echo ""
echo "${VOLATILEPATH}/ROOTfiles/ - ROOTfiles"
echo "${VOLATILEPATH}/OUTPUT/ - OUTPUT"
echo "${VOLATILEPATH}/REPORT_OUTPUT/ - REPORT_OUTPUT"
echo "${VOLATILEPATH}/HISTOGRAMS/ - HISTOGRAMS"
echo""
read -p "Proceed with sym link setup setup? <Y/N> " prompt2
if [[ $prompt2 == "y" || $prompt2 == "Y" || $prompt2 == "yes" || $prompt2 == "Yes" ]]; then
    echo "Creating sym links"
else 
    echo "Please update script with desired changes to sym links and re-run if desired"
    exit 1
fi

# Each loop checks if the link exists, if it doesn't, make it
# If it DOES, check it's not broken, if broken, replace, if not, just print that it exists
if [ ! -L "${UTILPATH}/ROOTfiles" ]; then
    ln -s "${VOLATILEPATH}/ROOTfiles/" "${UTILPATH}/ROOTfiles"
elif [ -L "${UTILPATH}/ROOTfiles" ]; then
    if [ ! -e "${UTILPATH}/ROOTfiles" ]; then
	echo "ROOTfiles sym link exits but is broken, replacing"
	rm "${UTILPATH}/ROOTfiles"
	ln -s "${VOLATILEPATH}/ROOTfiles/" "${UTILPATH}/ROOTfiles"
    else echo "${UTILPATH}/ROOTfiles sym link already exists and not broken"
    fi
fi

if [ ! -L "${UTILPATH}/OUTPUT" ]; then
    ln -s "${VOLATILEPATH}/OUTPUT/" "${UTILPATH}/OUTPUT"
elif [ -L "${UTILPATH}/OUTPUT" ]; then
    if [ ! -e "${UTILPATH}/OUTPUT" ]; then
	echo "OUTPUT sym link exits but is broken, replacing"
	rm "${UTILPATH}/OUTPUT"
	ln -s "${VOLATILEPATH}/OUTPUT/" "${UTILPATH}/OUTPUT"
    else echo "${UTILPATH}/OUTPUT sym link already exists and not broken"
    fi
fi

if [ ! -L "${UTILPATH}/REPORT_OUTPUT" ]; then
    ln -s "${VOLATILEPATH}/REPORT_OUTPUT/" "${UTILPATH}/REPORT_OUTPUT"
elif [ -L "${UTILPATH}/REPORT_OUTPUT" ]; then
    if [ ! -e "${UTILPATH}/REPORT_OUTPUT" ]; then
	echo "REPORT_OUTPUT sym link exits but is broken, replacing"
	rm "${UTILPATH}/REPORT_OUTPUT"
	ln -s "${VOLATILEPATH}/REPORT_OUTPUT/" "${UTILPATH}/REPORT_OUTPUT"
    else echo "${UTILPATH}/REPORT_OUTPUT sym link already exists and not broken"
    fi
fi

if [ ! -L "${UTILPATH}/HISTOGRAMS" ]; then
    ln -s "${VOLATILEPATH}/HISTOGRAMS/" "${UTILPATH}/HISTOGRAMS"
elif [ -L "${UTILPATH}/HISTOGRAMS" ]; then
    if [ ! -e "${UTILPATH}/HISTOGRAMS" ]; then
	echo "HISTOGRAMS sym link exits but is broken, replacing"
	rm "${UTILPATH}/HISTOGRAMS"
	ln -s "${VOLATILEPATH}/HISTOGRAMS/" "${UTILPATH}/HISTOGRAMS"
    else echo "${UTILPATH}/HISTOGRAMS sym link already exists and not broken"
    fi
fi

echo ""
echo "Directories and sym links for UTIL_PION now setup"

exit 0
