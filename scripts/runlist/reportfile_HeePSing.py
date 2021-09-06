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

Current=0
PS1=0
PS4=0
PS5=0
HMS_Rate=0
SHMS_Rate=0
Coin_Rate=0
Charge=0
Raw_HMS=0
Raw_SHMS=0
Raw_Coin=0
EDTM=0
Elec_Track=0

TestVar = 0 # Counter to check the right number of variables have been set, should get 10 items
for line in ReportFile:
    if "SW_BCM4A_Current" in line :
        Current = float(((line.split(":")[1]).strip()).split(" ")[0]) # Need to split on : delimiter to get number, then space to remove unit
        TestVar+=1
    if "SW_Ps1_factor" in line :
        PS1 = int((line.split(":"))[1])
        TestVar+=1
    if "SW_Ps3_factor" in line : # In the HeePSingles runs, we actually only care about EL-Reals so we pick up the PS3 rate, to keep the number of columns down, this is printed as "PS4"
        PS4 = int((line.split(":"))[1])
        TestVar+=1
    if "SW_Ps5_factor" in line :
        PS5 = int((line.split(":"))[1])
        TestVar+=1
    if "SW_HMS_EL-REAL_Trigger_Rate" in line :
        HMS_Rate = float(((line.split(":")[1]).strip()).split(" ")[0]) 
        TestVar+=1
    if "SW_SHMS_EL-REAL_Trigger_Rate" in line :
        SHMS_Rate = float(((line.split(":")[1]).strip()).split(" ")[0]) 
        TestVar+=1
    if "SW_COIN_Trigger_Rate" in line :
        COIN_Rate = float(((line.split(":")[1]).strip()).split(" ")[0]) 
        TestVar+=1
    if "SW_BCM4A_Beam_Cut_Charge" in line :
        Charge = float(((line.split(":")[1]).strip()).split(" ")[0]) 
        TestVar+=1
    if "SW_Accepted_HMS_Triggers" in line :
        Raw_HMS = float(((line.split(":")[1]).strip()).split(" ")[0]) 
        TestVar+=1
    if "SW_Accepted_SHMS_Triggers" in line :
        Raw_SHMS = float(((line.split(":")[1]).strip()).split(" ")[0]) 
        TestVar+=1
    if "SW_Accepted_COIN_Triggers" in line :
        Raw_Coin = float(((line.split(":")[1]).strip()).split(" ")[0]) 
        TestVar+=1
    if "SW_EDTM_Accepted_Triggers" in line :
        EDTM = float(((line.split(":")[1]).strip()).split(" ")[0]) 
        TestVar+=1
    if "SW_SHMS_Electron_Singles_TRACK_EFF" in line :
        Elec_Track = float(((line.split(":")[1]).strip()).split(" ")[0])  
        TestVar+=1

#Round these values to nearest 1000
Raw_HMS=int(round(Raw_HMS,-3)/1000)
Raw_SHMS=int(round(Raw_SHMS,-3)/1000)
Raw_Coin=int(round(Raw_Coin,-3)/1000)
EDTM=int(round(EDTM,-3)/1000)
        
if TestVar != 13 and TestVar > 13 :
    print(" !!! WARNING IN reportfile.py !!! \n More than expected matching entries found, some information may have been overwritten \n !!! WARNING IN reportfile.py !!!")
    RunListEntry=("%.3f,%i,%i,%i,%.3f,%.3f,%.3f,%.3f,%ik,%ik,%ik,%ik,%.3f" % (Current, PS1, PS4, PS5, HMS_Rate, SHMS_Rate, COIN_Rate, Charge, Raw_HMS, Raw_SHMS, Raw_Coin, EDTM, Elec_Track) )
elif TestVar != 13 and TestVar < 13 :
    print(" !!! WARNING IN reportfile.py !!! \n Less than expected matching entries found, some information may have not have been gathered \n !!! WARNING IN reportfile.py !!!")
    RunListEntry=("%.3f,%i,%i,%i,%.3f,%.3f,%.3f,%.3f,%ik,%ik,%ik,%ik,%.3f" % (Current, PS1, PS4, PS5, HMS_Rate, SHMS_Rate, COIN_Rate, Charge, Raw_HMS, Raw_SHMS, Raw_Coin, EDTM, Elec_Track) )
else :
    RunListEntry=("%.3f,%i,%i,%i,%.3f,%.3f,%.3f,%.3f,%ik,%ik,%ik,%ik,%.3f" % (Current, PS1, PS4, PS5, HMS_Rate, SHMS_Rate, COIN_Rate, Charge, Raw_HMS, Raw_SHMS, Raw_Coin, EDTM, Elec_Track) )
print(RunListEntry)

ReportFile.close()
