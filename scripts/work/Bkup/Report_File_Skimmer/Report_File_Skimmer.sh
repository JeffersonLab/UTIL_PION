#! /bin/bash

### Author info
### Script to skim over list of runs and grab info from their report files

RunList=$1
if [[ -z "$1" ]]; then
    echo "I need a list of run numbers to process!"
    echo "Please provide a run list as an input variable"
    exit 1
fi
OutputFile=$2
if [[ -z "$2" ]]; then
    echo "I need an output file name to write to!"
    echo "Please provide an output file name!"
    exit 2
fi
# Report file path, change as needed
ReportPath="/volatile/hallc/c-kaonlt/aliusman/REPORT_OUTPUT/Pass2_Files/"

### You should define your header and print it to the top of your new output file

# Loop over all report files for run numbers listed in input
while IFS='' read -r line || [[ -n "$line" ]]; do
    RunNum=$line
    ReportFile="${ReportPath}/replay_coin_production_${RunNum}_-1.report"
    ### Add a line here which checks the file exists, skip it if it doesn't exist. You probably also want to add it to a list to check later
    ### If the file exists, grab the variables you want with sed commands, set the output of the sed command to be the variable
    ### Once you have all of the variables you need, append it to the output file with an echo command via >>
    ### After the echo command, this report file should be done
done < "$RunList"
