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
COIN_Rate=0
Charge=0
Raw_HMS=0
Raw_SHMS=0
Raw_COIN=0
EDTM=0
Had_Track=0

TestVar = 0 # Counter to check the right number of variables have been set, should get 10 items
for line in ReportFile:
    if "SW_BCM4A_Current" in line :
        Current = float(((line.split(":")[1]).strip()).split(" ")[0]) # Need to split on : delimiter to get number, then space to remove unit
        TestVar+=1
        #print('Current', TestVar, "\n")
    if "SW_Ps1_factor" in line :
        PS1 = int((line.split(":"))[1])
        TestVar+=1
        #print('PS1', TestVar, "\n")
    if "SW_Ps4_factor" in line :
        PS4 = int((line.split(":"))[1])
        TestVar+=1
        #print('PS4', TestVar, "\n")
    if "SW_Ps5_factor" in line :
        PS5 = int((line.split(":"))[1])
        TestVar+=1
        #print('PS5', TestVar, "\n")
    if "SW_HMS_EL-REAL_Trigger_Rate" in line :
        HMS_Rate = float(((line.split(":")[1]).strip()).split(" ")[0]) 
        TestVar+=1
        #print('HMS_Rate', TestVar, "\n")
    if "SW_SHMS_3/4_Trigger_Rate" in line :
        SHMS_Rate = float(((line.split(":")[1]).strip()).split(" ")[0]) 
        TestVar+=1
        #print('SHMS_Rate', TestVar, "\n")
    if "SW_COIN_Trigger_Rate" in line :
        COIN_Rate = float(((line.split(":")[1]).strip()).split(" ")[0]) 
        TestVar+=1
        #print('COIN_Rate', TestVar, "\n")
    if "SW_BCM4A_Beam_Cut_Charge" in line :
        Charge = float(((line.split(":")[1]).strip()).split(" ")[0]) 
        TestVar+=1
        #print('Charge', TestVar, "\n")
    if "SW_Pre-Scaled_HMS_EL-REAL_Triggers" in line :
        Raw_HMS = float(((line.split(":")[1]).strip()).split(" ")[0]) 
        TestVar+=1
        #print('Raw_HMS', TestVar, "\n")
    if "SW_Pre-Scaled_SHMS_3/4_Triggers" in line :
        Raw_SHMS = float(((line.split(":")[1]).strip()).split(" ")[0]) 
        TestVar+=1
        #print('Raw_SHMS', TestVar, "\n")
    if "SW_Pre-Scaled_COIN_Triggers" in line :
        Raw_COIN = float(((line.split(":")[1]).strip()).split(" ")[0]) 
        TestVar+=1
        #print('Raw_COIN', TestVar, "\n")
    if "SW_EDTM_Accepted_Triggers" in line :
        EDTM = float(((line.split(":")[1]).strip()).split(" ")[0]) 
        TestVar+=1
        #print('EDTM', TestVar, "\n")
    if "SW_SHMS_Hadron_Coin_TRACK_EFF" in line :
        Had_Track = float(((line.split(":")[1]).strip()).split(" ")[0])  
        TestVar+=1
        #print('Had_Track', TestVar, "\n")

#Round these values to nearest 1000
Raw_HMS=abs(int(round(Raw_HMS,-3)/1000))
Raw_SHMS=abs(int(round(Raw_SHMS,-3)/1000))
Raw_COIN=abs(int(round(Raw_COIN,-3)/1000))
EDTM=int(round(EDTM,-3)/1000)

# Need to convert PS actual value which is read in above as PS1/PS4/PS5 to the PS value (-1 to 16) that is actually filled in the run list
psActual = [-1,1,2,3,5,9,17,33,65,129,257,513,1025,2049,4097,8193,16385,32769]
psValue = [-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]

for i,val in enumerate(psActual):
    if val == PS1:
        if(val == -1):
            PS1 = 0
        else:
            PS1 = int(psValue[i])
    if val == PS4:
        if(val == -1):
            PS4 = 0
        else:
            PS4 = int(psValue[i])
    if val == PS5:
        # if(val == -1):
        #     PS5 = 0
        # else:
        PS5 = int(psValue[i])

if TestVar != 13 and TestVar > 13 :
    print(" !!! WARNING IN reportfile.py !!! \n More than expected matching entries found, some information may have been overwritten \n !!! WARNING IN reportfile.py !!!")
    RunListEntry=("%.3f,%i,%i,%i,%.3f,%.3f,%.3f,%.3f,%ik,%ik,%ik,%ik,%.3f" % (Current, PS1, PS4, PS5, HMS_Rate, SHMS_Rate, COIN_Rate, Charge, Raw_HMS, Raw_SHMS, Raw_COIN, EDTM, Had_Track) )
elif TestVar != 13 and TestVar < 13 :
    print(" !!! WARNING IN reportfile.py !!! \n Less than expected matching entries found, some information may have not have been gathered \n !!! WARNING IN reportfile.py !!!")
    RunListEntry=("%.3f,%i,%i,%i,%.3f,%.3f,%.3f,%.3f,%ik,%ik,%ik,%ik,%.3f" % (Current, PS1, PS4, PS5, HMS_Rate, SHMS_Rate, COIN_Rate, Charge, Raw_HMS, Raw_SHMS, Raw_COIN, EDTM, Had_Track) )
else :
    RunListEntry=("%.3f,%i,%i,%i,%.3f,%.3f,%.3f,%.3f,%ik,%ik,%ik,%ik,%.3f" % (Current, PS1, PS4, PS5, HMS_Rate, SHMS_Rate, COIN_Rate, Charge, Raw_HMS, Raw_SHMS, Raw_COIN, EDTM, Had_Track) )
print(RunListEntry)

ReportFile.close()
