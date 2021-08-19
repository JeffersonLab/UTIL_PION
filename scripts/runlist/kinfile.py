#! /usr/bin/python

# Stephen JD Kay - 19/08/21 - University of Regina
# Script to grab info from standard.kinematics, for use in the run list shell script
# This script splits the standard.kinematics file into blocks, and loops over every line in a block. It strips out comments such that the first line is the run number range
# It finds the relevant run number range block, and gets info to feed into the run list
# It will send the info from the LAST matching block in the kinematics file

# Import relevant packages
import sys, math, os, subprocess

sys.path.insert(0, 'python/')

if len(sys.argv)-1!=2:
    print("!!!!! ERROR !!!!!\n Expected 2 arguments\n Usage is with - KinFilePath and run number \n!!!!! ERROR !!!!!")
    sys.exit(1)

KinFilePath = sys.argv[1]
RunNum = sys.argv[2]

KinFile = open(KinFilePath)
KinFileContent = KinFile.read()
KinFile.close()

SHMS_Angle = "ERROR"
SHMS_P = "ERROR"
HMS_Angle = "ERROR"
HMS_P = "ERROR"
EBeam = "ERROR"

TestVar = 0
# The loop here is explicitly written such that if there are multiple entrie, the values will be overwrriten
# The LAST matching block in the file is the one that will be used
for KinFileBlock in KinFileContent.split('\n\n'):
    nLines=0 # Counter for the number of lines in the block of text
    Lines =[]
    for KinFileLine in KinFileBlock.split('\n'):
        nLines+=1
        if not KinFileLine.startswith("#"): # If line does NOT start with a #, add it to our array
            Lines.append(KinFileLine)
    if nLines < 10: # If less than 10 lines, skip to next block
        continue
    # If it's an entry with a -, it's a range of run numbers, set the start and end accordingly    
    if "-" in Lines[0]:
        RunNumArr = Lines[0].split("-")
        RunStart = int(RunNumArr[0])
        RunEnd = int(RunNumArr[1])
    # If there's no -, it's a single line entry and run start and end are the same
    elif "-" not in Lines[0]:
        RunStart=int(Lines[0])
        RunEnd=int(Lines[0])
    # Check if the provided run number is in the run number range for the block, if it is, set the values
    if int(RunNum) in range (RunStart, RunEnd) :
        TestVar +=1
        for entry in Lines :
            if "ptheta_lab" in entry :
                SHMS_Angle = float((entry.split("="))[1])
            if "ppcentral" in entry :
                SHMS_P = float((entry.split("="))[1])
            if "htheta_lab" in entry :
                HMS_Angle = float((entry.split("="))[1])
            if "hpcentral" in entry :
                HMS_P = float((entry.split("="))[1])
            if "gpbeam" in entry :
                EBeam = float((entry.split("="))[1])
RunListEntry=("%2.3f,%2.3f,%2.3f,%2.3f,%2.3f" % (SHMS_Angle, SHMS_P, HMS_Angle, HMS_P, EBeam))
print(RunListEntry)
            
            
