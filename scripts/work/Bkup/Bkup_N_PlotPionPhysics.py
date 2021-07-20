#! /usr/bin/python

# 21/05/21 - Muhammad Junaid, University of Regina, Canada

# Python version of the pion plotting script. Now utilises uproot to select event of each type and writes them to a root file.
# Python should allow for easier reading of databases storing diferent variables.

###################################################################################################################################################

# Import relevant packages
import uproot as up
import numpy as np
import root_numpy as rnp
import pandas as pd
import root_pandas as rpd
import ROOT
import scipy
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import sys, math, os, subprocess
import array
sys.path.insert(0, 'python/')

##################################################################################################################################################

# Check the number of arguments provided to the script
if len(sys.argv)-1!=3:
    print("!!!!! ERROR !!!!!\n Expected 3 arguments\n Usage is with - ROOTfilePrefix RunNumber MaxEvents \n!!!!! ERROR !!!!!")
    sys.exit(1)

##################################################################################################################################################

# Input params - run number and max number of events
runNum = sys.argv[1]
MaxEvent = sys.argv[2]
ROOTPrefix = sys.argv[3]

USER = subprocess.getstatusoutput("whoami") # Grab user info for file finding
HOST = subprocess.getstatusoutput("hostname")

if ("farm" in HOST[1]):
#    REPLAYPATH = "/group/c-pionlt/USERS/%s/hallc_replay_lt" % USER[1]
    REPLAYPATH = "/group/c-kaonlt/USERS/%s/hallc_replay_lt" % USER[1]

elif ("qcd" in HOST[1]):
#    REPLAYPATH = "/group/c-pionlt/USERS/%s/hallc_replay_lt" % USER[1]
    REPLAYPATH = "/group/c-kaonlt/USERS/%s/hallc_replay_lt" % USER[1]

elif ("phys.uregina" in HOST[1]):
    REPLAYPATH = "/home/%s/work/JLab/hallc_replay_lt" % USER[1]

elif("skynet" in HOST[1]):
    REPLAYPATH = "/home/%s/Work/JLab/hallc_replay_lt" % USER[1]

#################################################################################################################################################

# Add more path setting as needed in a similar manner                                                                                                                                                          
OUTPATH = "%s/UTIL_PION/scripts/demo/OUTPUT/Analysis/PionLT" % REPLAYPATH        # Output folder location                                                                                                     
sys.path.insert(0, '%s/UTIL_PION/bin/python/' % REPLAYPATH)
import kaonlt as klt # Import kaonlt module, need the path setting line above prior to importing this                                                                                                         
print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER[1], HOST[1], REPLAYPATH))

#################################################################################################################################################

# Construct the name of the rootfile based upon the info we provided
rootName = "%s/UTIL_PION//scripts/demo/OUTPUT/Analysis/PionLT/%s_%s_%s.root" % (REPLAYPATH, runNum, MaxEvent, ROOTPrefix)     # Input file location and variables taking
print ("Attempting to process %s" %(rootName))
if os.path.exists(OUTPATH):
    if os.path.islink(OUTPATH):
        pass
    elif os.path.isdir(OUTPATH):
        pass
    else:
        print ("%s exists but is not a directory or sym link, check your directory/link and try again" % (OUTPATH))
        sys.exit(2)
else:
    print("Output path not found, please make a sym link or directory called OUTPUT in UTIL_PION/scripts/demo to store output")
    sys.exit(3)
print ("Attempting to process %s" %(rootName))
if os.path.isfile(rootName):
    print ("%s exists, attempting to process" % (rootName))
else:
    print ("%s not found - do you have the correct sym link/folder set up?" % (rootName))
    sys.exit(4)
print("Output path checks out, outputting to %s" % (OUTPATH))

###############################################################################################################################################

# Read stuff from the main event tree

infile = ROOT.TFile.Open(rootName, "READ")
Uncut_Pion_Events_tree = infile.Get("Uncut_Pion_Events")

TestBranch=Uncut_Pion_Events_tree.GetBranch('H_gtr_xp')

nEntries = Uncut_Pion_Events_tree.GetEntries()
for entryNum in range (0, nEntries):
    Uncut_Pion_Events_tree.GetEntry(entryNum)
    if entryNum > 100:
        break

#print(Uncut_Pion_Events_tree.GetEntries())
#print(TestBranch.GetEntries())
#for event in Uncut_Pion_Events_tree:
 #   print((Uncut_Pion_Events_tree.GetBranch('H_gtr_xp')).GetEntry(i))
  #  if i > 100:
   #     break

#Uncut_Pion_Events_tree = up.open(rootName)["Uncut_Pion_Events"]

#nEntries = Uncut_Pion_Events_tree.GetEntries()


# Variables from Uncut_Pion_Events TTree
#H_gtr_beta_pions_nocut = Uncut_Pion_Events_tree.array("H_gtr_beta")

# Branch
#t_tree.Branch("H_gtr_beta", H_gtr_beta)

###############################################################################################################################################

# Defining Histograms
H_gtr_beta_pions_Uncut = ROOT.TH1D("h1_H_gtr_beta_pions_Uncut", "H_gtr_beta", 220, 0.6, 1.4)

# Filling Histograms
#for entryNum in range(0, nEntries):
#	Uncut_Pion_Events_tree.GetEntry(entryNum)
#	H_gtr_beta_pions_Uncut.Fill(H_gtr_beta_pions_nocut)

#H_gtr_beta_pions_Uncut.Fill(H_gtr_beta_pions_nocut)

##############################################################################################################################################

#infile.Close()

#outHistFile = ROOT.TFile("outFileName.root", "RECREATE")
#outHistFile.cd()

# Writing Histograms
#H_gtr_beta_pions_uncut.Write()

#outHistFile.Close()

#############################################################################################################################################

def coin_pions():

    # Define the array of arrays containing the relevant HMS and SHMS info
    NoCut_COIN_Pions = [H_gtr_beta_pions_uncut]
    Uncut_COIN_Pions = [(H_gtr_beta_pions_uncut) for (H_gtr_beta_pions_uncut) in zip(*NoCut_COIN_Pions)]
    COIN_Pions = {
        "Uncut_Pion_Events" : Uncut_COIN_Pions,
        }

    return COIN_Pions

############################################################################################################################################

#def main():
#    COIN_Pion_Data = coin_pions()

#    COIN_Pion_Data_Header = ["H_gtr_beta_pions_uncut"]
    # Need to create a dict for all the branches we grab
#    data = {}

#    for d in (COIN_Pion_Data): # Convert individual dictionaries into a "dict of dicts"
 #       data.update(d)
 #       data_keys = list(data.keys()) # Create a list of all the keys in all dicts added above, each is an array of data

  #  for i in range (0, len(data_keys)):
   #     if("Pion" in data_keys[i]):
    #        DFHeader=list(COIN_Pion_Data_Header)
     #   else:
      #      continue
       #     # Uncomment the line below if you want .csv file output, WARNING the files can be very large and take a long time to process!
            #pd.DataFrame(data.get(data_keys[i])).to_csv("%s/%s_%s.csv" % (OUTPATH, data_keys[i], runNum), header=DFHeader, index=False) # Convert array to panda dataframe and write to csv with correct header
       # if (i == 0):
        #    pd.DataFrame(data.get(data_keys[i]), columns = DFHeader, index = None).to_root("%s/%s_%s_Output_Data.root" % (OUTPATH, runNum, MaxEvent), key ="%s" % data_keys[i])
        #elif (i != 0):
         #   pd.DataFrame(data.get(data_keys[i]), columns = DFHeader, index = None).to_root("%s/%s_%s_Output_Data.root" % (OUTPATH, runNum, MaxEvent), key ="%s" % data_keys[i], mode ='a')

#if __name__ == '__main__':
#    main()
#print ("Processing Complete")
