#! /usr/bin/python

# 09/10/20 - Stephen Kay, University of Regina

# A short python demo script demonstrating opening up a root file, using uproot to grab some info and then save it as a new rootfile

# Import relevant packages, there's more than we need here but there's no harm in that
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

sys.path.insert(0, 'python/')
# Check the number of arguments provided to the script
if len(sys.argv)-1!=3:
    print("!!!!! ERROR !!!!!\n Expected 3 arguments\n Usage is with - ROOTfilePrefix RunNumber MaxEvents \n!!!!! ERROR !!!!!")
    sys.exit(1)
# Input params - run number and max number of events
ROOTPrefix = sys.argv[1]
runNum = sys.argv[2]
MaxEvent = sys.argv[3]
print (ROOTPrefix, runNum, MaxEvent)

USER = subprocess.getstatusoutput("whoami") # Grab user info for file finding
HOST = subprocess.getstatusoutput("hostname")
if ("farm" in HOST[1]):
    REPLAYPATH = "/group/c-kaonlt/USERS/%s/hallc_replay_lt" % USER[1]
elif ("qcd" in HOST[1]):
    REPLAYPATH = "/group/c-kaonlt/USERS/%s/hallc_replay_lt" % USER[1]
elif ("phys.uregina" in HOST[1]):
    REPLAYPATH = "/home/%s/work/JLab/hallc_replay_lt" % USER[1]
elif ("skynet" in HOST[1]):
    REPLAYPATH = "/home/%s/Work/JLab/hallc_replay_lt" % USER[1]

# Add more path setting as needed in a similar manner
OUTPATH = "%s/UTIL_KAONLT/scripts/demo/OUTPUT" % REPLAYPATH
CUTPATH = "%s/UTIL_KAONLT/DB/CUTS" % REPLAYPATH
sys.path.insert(0, '%s/UTIL_KAONLT/bin/python/' % REPLAYPATH)
import kaonlt as klt # Import kaonlt module, need the path setting line above prior to importing this

print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER[1], HOST[1], REPLAYPATH))
# Construct the name of the rootfile based upon the info we provided
rootName = "%s/UTIL_KAONLT/ROOTfiles/%s_%s_%s.root" % (REPLAYPATH, ROOTPrefix, runNum, MaxEvent)

# Read stuff from the main event tree, here we're just going to get some quantities for the acceptance for the HMS/SHMS
e_tree = up.open(rootName)["T"]
# HMS info
H_gtr_beta = e_tree.array("H.gtr.beta")
H_gtr_xp = e_tree.array("H.gtr.th") # xpfp -> Theta
H_gtr_yp = e_tree.array("H.gtr.ph") # ypfp -> Phi
H_gtr_dp = e_tree.array("H.gtr.dp")
# SHMS info
P_gtr_beta = e_tree.array("P.gtr.beta")
P_gtr_xp = e_tree.array("P.gtr.th") # xpfp -> Theta
P_gtr_yp = e_tree.array("P.gtr.ph") # ypfp -> Phi
P_gtr_p = e_tree.array("P.gtr.p")
P_gtr_dp = e_tree.array("P.gtr.dp")

# Define a function to return a dictionary of the events we want
# Arrays we generate in our dict should all be of the same length (in terms of # elements in the array) to keep things simple
def All_events(): 
    # Define the array of arrays containing the relevant HMS and SHMS info
    NoCut_Events = [H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp] # Create a LIST of the arrays we want
    # Turn our list explicitly into an NP array, the names here don't matter too much so long as they're the same on each side - they match up to the 1st/2nd...nth element of our list basically
    Events_Info = [(HBeta, Hxp, Hyp, Hdel, PBeta, Pxp, Pyp, PP, PDel) for (HBeta, Hxp, Hyp, Hdel, PBeta, Pxp, Pyp, PP, PDel) in zip(*NoCut_Events)]
    
    All_Events = { # Create a dictionary containg our array
        # Name of the element in our dictionary - this is important!
        # The elements (dict keys) we enter here are the "trees" in the root file we make at the end, the branches of the tree are the elements of the array above (which are themselves, arrays)
        # Our produced root file will have trees named according to what we enter here
        "All_Events" : Events_Info,
    }

    return All_Events # Return a dictionary

# The name here is a little misleading, there are no cuts so "HMS_Events" here just refers to the fact that we only have HMS info in this dict
def HMS_events(): 
    NoCut_Events_HMS = [H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp]
    HMS_Events_Info = [(HBeta, Hxp, Hyp, Hdel) for (HBeta, Hxp, Hyp, Hdel) in zip(*NoCut_Events_HMS)]

    HMS_Events = {
        "HMS_Events" : HMS_Events_Info,
    }

    return HMS_Events

# See comment above on the "HMS events"
def SHMS_events(): 
    NoCut_Events_SHMS = [P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp]
    SHMS_Events_Info = [(PBeta, Pxp, Pyp, PP, PDel) for (PBeta, Pxp, Pyp, PP, PDel) in zip(*NoCut_Events_SHMS)]

    SHMS_Events = {
        "SHMS_Events" : SHMS_Events_Info,
    }

    return SHMS_Events

def main():
    # Run our functions and get a dict from each
    All_Events_Data = All_events()    
    HMS_Events_Data = HMS_events()
    SHMS_Events_Data = SHMS_events()
    
    # This is just the list of branches we use from the initial root file for each dict
    # They're the "headers" of the data frame we create - i.e. they're going to be the names of the tree branches in our new root file
    # Note - I don't like re-defining this here as it's very prone to errors if you included (or removed something) earlier but didn't modify it here
    All_Data_Header = ["H_gtr_beta","H_gtr_xp","H_gtr_yp","H_gtr_dp","P_gtr_beta","P_gtr_xp","P_gtr_yp","P_gtr_p","P_gtr_dp"] 
    HMS_Data_Header = ["H_gtr_beta","H_gtr_xp","H_gtr_yp","H_gtr_dp"]
    SHMS_Data_Header = ["P_gtr_beta","P_gtr_xp","P_gtr_yp","P_gtr_p","P_gtr_dp"]
    data = {} # Create an empty dictionary

    for d in (All_Events_Data, HMS_Events_Data, SHMS_Events_Data): # Convert individual dictionaries into a "dict of dicts"
        data.update(d) # For every dictionary we give above, add its keys to the new dict
        data_keys = list(data.keys()) # Create a list of all the keys in all dicts added above, each is an array of data

    for i in range (0, len(data_keys)):
        # Set the headers for our data frame
        if ("All_" in data_keys[i]):
            DFHeader=list(All_Data_Header)
        elif ("HMS_" in data_keys[i]):
            if("SHMS_" in data_keys[i]):
                DFHeader=list(SHMS_Data_Header)
            else:
                DFHeader=list(HMS_Data_Header)
        else:
            continue
        if (i == 0): # For the first case, start writing to file
            pd.DataFrame(data.get(data_keys[i]), columns = DFHeader, index = None).to_root("%s/%s_%s_Demo1_Data.root" % (OUTPATH, runNum, MaxEvent), key ="%s" % data_keys[i])
        elif (i != 0): # For any but the first case, append it to our file
            pd.DataFrame(data.get(data_keys[i]), columns = DFHeader, index = None).to_root("%s/%s_%s_Demo1_Data.root" % (OUTPATH, runNum, MaxEvent), key ="%s" % data_keys[i], mode ='a') 
                    
if __name__ == '__main__': # Execute our "main" function that we defined above
    main()
