#! /usr/bin/python

# Stephen JD Kay - 19/08/21 - University of Regina
# Script to grab info from a report file, for use in the run list shell script
# This script scans every line of the report file to find the correct info

# Import relevant packages
import sys, math, os, subprocess

sys.path.insert(0, 'python/')

if len(sys.argv)-1!=1:
    print("!!!!! ERROR !!!!!\n Expected 1 argument\n Usage is with - ReportFilePath \n!!!!! ERROR !!!!!")
    sys.exit(1)

ReportFilePath = sys.argv[1]

ReportFile = open(ReportFilePath)

Current="ERROR"
PS1="ERROR"
PS4="ERROR"
PS5="ERROR"
HMS_Rate="ERROR"
SHMS_Rate="ERROR"
Coin_Rate="ERROR"
Charge="ERROR"
Raw_Coin="ERROR"
Had_Track="ERROR"

TestVar = 0 # Counter to check the right number of variables have been set, should get 10 items
for line in ReportFile:
    if "BCM4A Current" in line :
        Current = float(((line.split(":")[1]).strip()).split(" ")[0]) # Need to split on : delimiter to get number, then space to remove unit
        TestVar+=1
    if "Ps1_factor" in line :
        PS1 = int((line.split(":"))[1])
        TestVar+=1
    if "Ps4_factor" in line :
        PS4 = int((line.split(":"))[1])
        TestVar+=1
    if "Ps6_factor" in line :
        PS6 = int((line.split(":"))[1])
        TestVar+=1
    if "HMS EL-REAL Trigger Rate" in line :
        HMS_Rate = float(((line.split(":")[1]).strip()).split(" ")[0]) 
        TestVar+=1
    if "SHMS 3/4 Trigger Rate" in line :
        SHMS_Rate = float(((line.split(":")[1]).strip()).split(" ")[0]) 
        TestVar+=1
    if "COIN Trigger Rate" in line :
        COIN_Rate = float(((line.split(":")[1]).strip()).split(" ")[0]) 
        TestVar+=1
    if "BCM4A Beam Cut Charge" in line :
        Charge = float(((line.split(":")[1]).strip()).split(" ")[0]) 
        TestVar+=1
    if "Accepted COIN Triggers" in line :
        Raw_Coin = float(((line.split(":")[1]).strip()).split(" ")[0]) 
        TestVar+=1
    if "SHMS Hadron Singles TRACK EFF" in line :
        Had_Track = float(((line.split(":")[1]).strip()).split(" ")[0])  
        TestVar+=1

if TestVar != 10 :
    print(" !!! WARNING IN reportfile.py !!! \n More than expected matching entries found, some information may have been overwritten \n !!! WARNING IN reportfile.py !!!")

# # The loop here is explicitly written such that if there are multiple entrie, the values will be overwrriten
# # The LAST matching block in the file is the one that will be used
# for KinFileBlock in KinFileContent.split('\n\n'):
#     nLines=0 # Counter for the number of lines in the block of text
#     Lines =[]
#     for KinFileLine in KinFileBlock.split('\n'):
#         nLines+=1
#         if not KinFileLine.startswith("#"): # If line does NOT start with a #, add it to our array
#             Lines.append(KinFileLine)
#     if nLines < 10: # If less than 10 lines, skip to next block
#         continue
#     # If it's an entry with a -, it's a range of run numbers, set the start and end accordingly    
#     if "-" in Lines[0]:
#         RunNumArr = Lines[0].split("-")
#         RunStart = int(RunNumArr[0])
#         RunEnd = int(RunNumArr[1])
#     # If there's no -, it's a single line entry and run start and end are the same
#     elif "-" not in Lines[0]:
#         RunStart=int(Lines[0])
#         RunEnd=int(Lines[0])
#     # Check if the provided run number is in the run number range for the block, if it is, set the values
#     if int(RunNum) in range (RunStart, RunEnd) :
#         TestVar +=1
#         for entry in Lines :
#             if "ptheta_lab" in entry :
#                 SHMS_Angle = float((entry.split("="))[1])
#             if "ppcentral" in entry :
#                 SHMS_P = float((entry.split("="))[1])
#             if "htheta_lab" in entry :
#                 HMS_Angle = float((entry.split("="))[1])
#             if "hpcentral" in entry :
#                 HMS_P = float((entry.split("="))[1])
#             if "gpbeam" in entry :
#                 EBeam = float((entry.split("="))[1])
# RunListEntry=("%2.3f,%2.3f,%2.3f,%2.3f,%2.3f" % (SHMS_Angle, SHMS_P, HMS_Angle, HMS_P, EBeam))
# print(RunListEntry)
ReportFile.close()
