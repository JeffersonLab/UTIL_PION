# /usr/bin/python

# 09/10/20 - Stephen Kay, University of Regina

# A short python demo script demonstrating opening up a root file, using uproot to grab some info and then save it as a new rootfile
# This time, we also define and apply some cuts to our SHMS events

#########################################################################################################################
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

sys.path.insert(0, 'python/')

##########################################################################################################################

# Check the number of arguments provided to the script
if len(sys.argv)-1!=3:
    print("!!!!! ERROR !!!!!\n Expected 3 arguments\n Usage is with - ROOTfilePrefix RunNumber MaxEvents \n!!!!! ERROR !!!!!")
    sys.exit(1)

##########################################################################################################################

# Input params - run number and max number of events
ROOTPrefix = sys.argv[1]
runNum = sys.argv[2]
MaxEvent = sys.argv[3]

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

###########################################################################################################################

# Add more path setting as needed in a similar manner
OUTPATH = "%s/UTIL_PION/scripts/demo/OUTPUT" % REPLAYPATH        # Output folder location
sys.path.insert(0, '%s/UTIL_PION/bin/python/' % REPLAYPATH)

import kaonlt as klt # Import kaonlt module, need the path setting line above prior to importing this

print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER[1], HOST[1], REPLAYPATH))

############################################################################################################################

# Construct the name of the rootfile based upon the info we provided
rootName = "%s/UTIL_PION/ROOTfiles/%s_%s_%s.root" % (REPLAYPATH, ROOTPrefix, runNum, MaxEvent)     # Input file location and variables taking
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

##############################################################################################################################

# Read stuff from the main event tree, here we're just going to get some quantities for the acceptance for the HMS/SHMS
e_tree = up.open(rootName)["T"]

# HMS info
#H_gtr_beta = e_tree.array("H.gtr.beta")     # Beta is velocity of particle between pairs of hodoscopes
H_gtr_xp = e_tree.array("H.gtr.th")         # xpfp -> Theta
H_gtr_yp = e_tree.array("H.gtr.ph")         # ypfp -> Phi
H_gtr_dp = e_tree.array("H.gtr.dp")         # dp is Delta

H_hod_goodscinhit = e_tree.array("H.hod.goodscinhit")       #
H_hod_goodstarttime = e_tree.array("H.hod.goodstarttime")   #
H_hod_betanotrack = e_tree.array("H.hod.betanotrack")       #
H_hgcer_npeSum = e_tree.array("H.cer.npeSum")             #
H_cal_etotnorm = e_tree.array("H.cal.etotnorm")             #

H_dc_ntrack = e_tree.array("H.dc.ntrack")

# SHMS info
#P_gtr_beta = e_tree.array("P.gtr.beta")
P_gtr_xp = e_tree.array("P.gtr.th") # xpfp -> Theta
P_gtr_yp = e_tree.array("P.gtr.ph") # ypfp -> Phi
#P_gtr_p = e_tree.array("P.gtr.p")   # p is momentum
P_gtr_dp = e_tree.array("P.gtr.dp")

P_hod_goodscinhit = e_tree.array("P.hod.goodscinhit")       #
P_hod_goodstarttime = e_tree.array("P.hod.goodstarttime")   #
P_hod_betanotrack = e_tree.array("P.hod.betanotrack")       #
P_hgcer_npeSum = e_tree.array("P.hgcer.npeSum")             #
P_aero_npeSum = e_tree.array("P.aero.npeSum")               #
P_cal_etotnorm = e_tree.array("P.cal.etotnorm")             #
P_dc_Ch1_nhit = e_tree.array("P.dc.Ch1.nhit")               #
P_dc_Ch2_nhit = e_tree.array("P.dc.Ch2.nhit")               #
P_dc_x_fp = e_tree.array("P.dc.x_fp")                       #
P_dc_y_fp = e_tree.array("P.dc.y_fp")                       #
P_dc_xp_fp = e_tree.array("P.dc.xp_fp")                     #
P_dc_yp_fp = e_tree.array("P.dc.yp_fp")                     #
XDip = P_dc_x_fp - P_dc_xp_fp*307
YDip = P_dc_y_fp - P_dc_yp_fp*307

P_dc_ntrack = e_tree.array("P.dc.ntrack")                   #

P_raw_ctime = e_tree.array("T.coin.pTRIG1_ROC2_tdcTimeRaw") #
H_raw_ctime = e_tree.array("T.coin.pTRIG4_ROC2_tdcTimeRaw") 
C_raw_ctime = P_raw_ctime - H_raw_ctime

r = klt.pyRoot()

# Specify the file which contains the cuts we want to use
#fout = '%s/UTIL_PION/DB/CUTS/run_type/demo.cuts' % REPLAYPATH

# read in cuts file and make dictionary
#c = klt.pyPlot(REPLAYPATH,DEBUG=True) # Switch False to True to enable DEBUG mode
#readDict = c.read_dict(fout,runNum)

################################################################################################################################

# This method calls several methods in kaonlt package. It is required to create properly formated
# dictionaries. The evaluation must be in the analysis script because the analysis variables (i.e. the
# leaves of interest) are not defined in the kaonlt package. This makes the system more flexible
# overall, but a bit more cumbersome in the analysis script. Perhaps one day a better solution will be
# implimented.
def make_cutDict(cut,inputDict=None):

    global c

    c = klt.pyPlot(REPLAYPATH,readDict)
    x = c.w_dict(cut)
    print("%s" % cut)
    print("x ", x)
    
    if inputDict == None:
        inputDict = {}
        
    for key,val in readDict.items():
        if key == cut:
            inputDict.update({key : {}})

    for i,val in enumerate(x):
        tmp = x[i]
        if tmp == "":
            continue
        else:
            inputDict[cut].update(eval(tmp))
        
    return inputDict

#################################################################################################################################
#Double_t XDip YDip
# [
#  XDip = P_dc_x_fp - P_dc_xp_fp*307,
#  YDip = P_dc_x_fp - P_dc_xp_fp*307,
# ]
# Add the cuts that we want to use from our specified file to the cut dictionary, note, we're only adding two of our three defined cuts to our cut dict
#cutDict = make_cutDict("Demo2Cut1")
#cutDict = make_cutDict("Demo2Cut2", cutDict)
#cutDict = make_cutDict("Demo2Cut3", cutDict)
#c = klt.pyPlot(REPLAYPATH,cutDict)

# Define a function to return a dictionary of the events we want
# Arrays we generate in our dict should all be of the same length (in terms of # elements in the array) to keep things simple
def All_events(): 
    # Define the array of arrays containing the relevant HMS and SHMS info

#    NoCut_Events = [H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp] # Create a LIST of the arrays we want
 
    NoCut_Events = [H_hod_goodscinhit, H_hod_goodstarttime, H_hod_betanotrack, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hgcer_npeSum, H_cal_etotnorm, H_dc_ntrack, P_hod_goodscinhit, P_hod_goodstarttime, P_hod_betanotrack, P_gtr_xp, P_gtr_yp, P_gtr_dp, P_hgcer_npeSum, P_aero_npeSum, P_cal_etotnorm, P_dc_Ch1_nhit, P_dc_Ch2_nhit, P_dc_x_fp, P_dc_y_fp, P_dc_xp_fp, P_dc_yp_fp, XDip, YDip, P_dc_ntrack] # Create a LIST of the arrays we want

   # Turn our list explicitly into an NP array, the names here don't matter too much so long as they're the same on each side - they match up to the 1st/2nd...nth element of our list basically
#    Events_Info = [(HBeta, Hxp, Hyp, Hdel, PBeta, Pxp, Pyp, PP, PDel) for (HBeta, Hxp, Hyp, Hdel, PBeta, Pxp, Pyp, PP, PDel) in zip(*NoCut_Events)]
 
    Events_Info = [(H_hod_goodscinhit, H_hod_goodstarttime, H_hod_betanotrack, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hgcer_npeSum, H_cal_etotnorm, H_dc_ntrack, P_hod_goodscinhit, P_hod_goodstarttime, P_hod_betanotrack, P_gtr_xp, P_gtr_yp, P_gtr_dp, P_hgcer_npeSum, P_aero_npeSum, P_cal_etotnorm, P_dc_Ch1_nhit, P_dc_Ch2_nhit, P_dc_x_fp, P_dc_y_fp, P_dc_xp_fp, P_dc_yp_fp, XDip, YDip, P_dc_ntrack) for (H_hod_goodscinhit, H_hod_goodstarttime, H_hod_betanotrack, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hgcer_npeSum, H_cal_etotnorm, H_dc_ntrack, P_hod_goodscinhit, P_hod_goodstarttime, P_hod_betanotrack, P_gtr_xp, P_gtr_yp, P_gtr_dp, P_hgcer_npeSum, P_aero_npeSum, P_cal_etotnorm, P_dc_Ch1_nhit, P_dc_Ch2_nhit, P_dc_x_fp, P_dc_y_fp, P_dc_xp_fp, P_dc_yp_fp, XDip, YDip, P_dc_ntrack) in zip(*NoCut_Events)]
   
    All_Events = { # Create a dictionary containg our array
        # Name of the element in our dictionary - this is important!
        # Our produced root file will have trees named according to what we enter here
        "All_Events" : Events_Info,
    }

    return All_Events # Return a dictionary

def SHMS_Only_events():

    NoCut_Events_SHMS_Only = [P_hod_goodscinhit, P_hod_goodstarttime, P_hod_betanotrack, P_gtr_xp, P_gtr_yp, P_gtr_dp, P_hgcer_npeSum, P_aero_npeSum, P_cal_etotnorm, P_dc_Ch1_nhit, P_dc_Ch2_nhit, P_dc_x_fp, P_dc_y_fp, P_dc_xp_fp, P_dc_yp_fp, XDip, YDip, P_dc_ntrack]
    SHMS_Only_Events_Info = [(P_hod_goodscinhit, P_hod_goodstarttime, P_hod_betanotrack, P_gtr_xp, P_gtr_yp, P_gtr_dp, P_hgcer_npeSum, P_aero_npeSum, P_cal_etotnorm, P_dc_Ch1_nhit, P_dc_Ch2_nhit, P_dc_x_fp, P_dc_y_fp, P_dc_xp_fp, P_dc_yp_fp, XDip, YDip, P_dc_ntrack) for (P_hod_goodscinhit, P_hod_goodstarttime, P_hod_betanotrack, P_gtr_xp, P_gtr_yp, P_gtr_dp, P_hgcer_npeSum, P_aero_npeSum, P_cal_etotnorm, P_dc_Ch1_nhit, P_dc_Ch2_nhit, P_dc_x_fp, P_dc_y_fp, P_dc_xp_fp, P_dc_yp_fp, XDip, YDip, P_dc_ntrack) in zip(*NoCut_Events_SHMS_Only)
    if P_hod_goodscinhit == 1
    if P_hod_goodstarttime == 1
    if P_hod_betanotrack > 0.4 and P_hod_betanotrack < 1.6
    if P_hgcer_npeSum > 1.5
    if P_aero_npeSum > 1.5
    if P_cal_etotnorm > 0.05 and P_cal_etotnorm < 0.6
    ]

    SHMS_Only_Events = {
        "SHMS_Only_Events" : SHMS_Only_Events_Info,
    }

    return SHMS_Only_Events


def Raw_Ctime_events():

    NoCut_Events_Raw_Ctime = [C_raw_ctime, P_hod_goodscinhit, H_hod_goodscinhit]
    Raw_Ctime_Events_Info = [(C_raw_ctime, P_hod_goodscinhit, H_hod_goodscinhit) for (C_raw_ctime, P_hod_goodscinhit, H_hod_goodscinhit) in zip(*NoCut_Events_Raw_Ctime)
    
    if P_hod_goodscinhit == 1
    if H_hod_goodscinhit == 1
    ]

    Raw_Ctime_Events = {
        "Raw_Ctime_Events" : Raw_Ctime_Events_Info,
    }

    return Raw_Ctime_Events

# The name here is a little misleading, there are no cuts so "HMS_Events" here just refers to the fact that we only have HMS info in this dict
def COIN_events(): 
#    NoCut_Events_HMS = [H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp]
#    HMS_Events_Info = [(HBeta, Hxp, Hyp, Hdel) for (HBeta, Hxp, Hyp, Hdel) in zip(*NoCut_Events_HMS)]

    NoCut_Events_COIN = [H_hod_goodscinhit, H_hod_goodstarttime, H_hod_betanotrack, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hgcer_npeSum, H_cal_etotnorm, H_dc_ntrack, P_hod_goodscinhit, P_hod_goodstarttime, P_hod_betanotrack, P_gtr_xp, P_gtr_yp, P_gtr_dp, P_hgcer_npeSum, P_aero_npeSum, P_cal_etotnorm, P_dc_Ch1_nhit, P_dc_Ch2_nhit, P_dc_x_fp, P_dc_y_fp, P_dc_xp_fp, P_dc_yp_fp, XDip, YDip, P_dc_ntrack]
    COIN_Events_Info = [(H_hod_goodscinhit, H_hod_goodstarttime, H_hod_betanotrack, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hgcer_npeSum, H_cal_etotnorm, H_dc_ntrack, P_hod_goodscinhit, P_hod_goodstarttime, P_hod_betanotrack, P_gtr_xp, P_gtr_yp, P_gtr_dp, P_hgcer_npeSum, P_aero_npeSum, P_cal_etotnorm, P_dc_Ch1_nhit, P_dc_Ch2_nhit, P_dc_x_fp, P_dc_y_fp, P_dc_xp_fp, P_dc_yp_fp, XDip, YDip, P_dc_ntrack) for (H_hod_goodscinhit, H_hod_goodstarttime, H_hod_betanotrack, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hgcer_npeSum, H_cal_etotnorm, H_dc_ntrack, P_hod_goodscinhit, P_hod_goodstarttime, P_hod_betanotrack, P_gtr_xp, P_gtr_yp, P_gtr_dp, P_hgcer_npeSum, P_aero_npeSum, P_cal_etotnorm, P_dc_Ch1_nhit, P_dc_Ch2_nhit, P_dc_x_fp, P_dc_y_fp, P_dc_xp_fp, P_dc_yp_fp, XDip, YDip, P_dc_ntrack) in zip(*NoCut_Events_COIN)

    if H_hod_goodscinhit == 1
    if H_hod_goodstarttime == 1
    if H_hod_betanotrack > 0.6 and H_hod_betanotrack < 1.4
    if H_hgcer_npeSum > 0.5
    if H_cal_etotnorm > 0.6
    if P_hod_goodscinhit == 1
    if P_hod_goodstarttime == 1
    if P_hod_betanotrack > 0.4 and P_hod_betanotrack < 1.6
    if P_hgcer_npeSum > 1.5
    if P_aero_npeSum > 1.5
    if P_cal_etotnorm > 0.05 and P_cal_etotnorm < 0.6
    ]

    COIN_Events = {
        "COIN_Events" : COIN_Events_Info,
    }

    return COIN_Events

def HMS_SHMS_Track1_events(): 

    NoCut_Events_HMS_SHMS_Track1 = [H_hod_goodscinhit, H_hod_goodstarttime, H_hod_betanotrack, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hgcer_npeSum, H_cal_etotnorm, H_dc_ntrack, P_hod_goodscinhit, P_hod_goodstarttime, P_hod_betanotrack, P_gtr_xp, P_gtr_yp, P_gtr_dp, P_hgcer_npeSum, P_aero_npeSum, P_cal_etotnorm, P_dc_Ch1_nhit, P_dc_Ch2_nhit, P_dc_x_fp, P_dc_y_fp, P_dc_xp_fp, P_dc_yp_fp, XDip, YDip, P_dc_ntrack]
    HMS_SHMS_Track1_Events_Info = [(H_hod_goodscinhit, H_hod_goodstarttime, H_hod_betanotrack, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hgcer_npeSum, H_cal_etotnorm, H_dc_ntrack, P_hod_goodscinhit, P_hod_goodstarttime, P_hod_betanotrack, P_gtr_xp, P_gtr_yp, P_gtr_dp, P_hgcer_npeSum, P_aero_npeSum, P_cal_etotnorm, P_dc_Ch1_nhit, P_dc_Ch2_nhit, P_dc_x_fp, P_dc_y_fp, P_dc_xp_fp, P_dc_yp_fp, XDip, YDip, P_dc_ntrack) for (H_hod_goodscinhit, H_hod_goodstarttime, H_hod_betanotrack, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hgcer_npeSum, H_cal_etotnorm, H_dc_ntrack, P_hod_goodscinhit, P_hod_goodstarttime, P_hod_betanotrack, P_gtr_xp, P_gtr_yp, P_gtr_dp, P_hgcer_npeSum, P_aero_npeSum, P_cal_etotnorm, P_dc_Ch1_nhit, P_dc_Ch2_nhit, P_dc_x_fp, P_dc_y_fp, P_dc_xp_fp, P_dc_yp_fp, XDip, YDip, P_dc_ntrack) in zip(*NoCut_Events_HMS_SHMS_Track1)

    if H_hod_goodscinhit == 1
    if H_hod_goodstarttime == 1
    if H_hod_betanotrack > 0.6 and H_hod_betanotrack < 1.4
    if H_hgcer_npeSum > 0.5
    if H_cal_etotnorm > 0.6
    if P_hod_goodscinhit == 1
    if P_hod_goodstarttime == 1
    if P_hod_betanotrack > 0.4 and P_hod_betanotrack < 1.6
    if P_hgcer_npeSum > 1.5
    if P_aero_npeSum > 1.5
    if P_cal_etotnorm > 0.05 and P_cal_etotnorm < 0.6
    if P_dc_ntrack == 1
    if abs(P_gtr_dp) > 50
    ]

    HMS_SHMS_Track1_Events = {
        "HMS_SHMS_Track1_Events" : HMS_SHMS_Track1_Events_Info,
    }

    return HMS_SHMS_Track1_Events
   
def HMS_SHMS_Trackn_events():

    NoCut_Events_HMS_SHMS_Trackn = [H_hod_goodscinhit, H_hod_goodstarttime, H_hod_betanotrack, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hgcer_npeSum, H_cal_etotnorm, H_dc_ntrack, P_hod_goodscinhit, P_hod_goodstarttime, P_hod_betanotrack, P_gtr_xp, P_gtr_yp, P_gtr_dp, P_hgcer_npeSum, P_aero_npeSum, P_cal_etotnorm, P_dc_Ch1_nhit, P_dc_Ch2_nhit, P_dc_x_fp, P_dc_y_fp, P_dc_xp_fp, P_dc_yp_fp, XDip, YDip, P_dc_ntrack]
    HMS_SHMS_Trackn_Events_Info = [(H_hod_goodscinhit, H_hod_goodstarttime, H_hod_betanotrack, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hgcer_npeSum, H_cal_etotnorm, H_dc_ntrack, P_hod_goodscinhit, P_hod_goodstarttime, P_hod_betanotrack, P_gtr_xp, P_gtr_yp, P_gtr_dp, P_hgcer_npeSum, P_aero_npeSum, P_cal_etotnorm, P_dc_Ch1_nhit, P_dc_Ch2_nhit, P_dc_x_fp, P_dc_y_fp, P_dc_xp_fp, P_dc_yp_fp, XDip, YDip, P_dc_ntrack) for (H_hod_goodscinhit, H_hod_goodstarttime, H_hod_betanotrack, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hgcer_npeSum, H_cal_etotnorm, H_dc_ntrack, P_hod_goodscinhit, P_hod_goodstarttime, P_hod_betanotrack, P_gtr_xp, P_gtr_yp, P_gtr_dp, P_hgcer_npeSum, P_aero_npeSum, P_cal_etotnorm, P_dc_Ch1_nhit, P_dc_Ch2_nhit, P_dc_x_fp, P_dc_y_fp, P_dc_xp_fp, P_dc_yp_fp, XDip, YDip, P_dc_ntrack) in zip(*NoCut_Events_HMS_SHMS_Trackn)

    if H_hod_goodscinhit == 1
    if H_hod_goodstarttime == 1
    if H_hod_betanotrack > 0.6 and H_hod_betanotrack < 1.4
    if H_hgcer_npeSum > 0.5
    if H_cal_etotnorm > 0.6
    if P_hod_goodscinhit == 1
    if P_hod_goodstarttime == 1
    if P_hod_betanotrack > 0.4 and P_hod_betanotrack < 1.6
    if P_hgcer_npeSum > 1.5
    if P_aero_npeSum > 1.5
    if P_cal_etotnorm > 0.05 and P_cal_etotnorm < 0.6
    if P_dc_ntrack > 1
    if abs(P_gtr_dp) > 50
    ]

    HMS_SHMS_Trackn_Events = {
        "HMS_SHMS_Trackn_Events" : HMS_SHMS_Trackn_Events_Info,
    }

    return HMS_SHMS_Trackn_Events


# See comment above on the "HMS events"

#####################################################################################################################################

def main():
    # Run our functions and get a dict from each
    All_Events_Data = All_events()    
    SHMS_Only_Events_Data = SHMS_Only_events()
    Raw_Ctime_Events_Data = Raw_Ctime_events()
    COIN_Events_Data = COIN_events()
    HMS_SHMS_Track1_Events_Data = HMS_SHMS_Track1_events()
    HMS_SHMS_Trackn_Events_Data = HMS_SHMS_Trackn_events()

    # This is just the list of branches we use from the initial root file for each dict
    # They're the "headers" of the data frame we create - i.e. they're going to be the branches in our new root file
    # Note - I don't like re-defining this here as it's very prone to errors if you included (or removed something) earlier but didn't modify it here
#    All_Data_Header = ["H_gtr_beta","H_gtr_xp","H_gtr_yp","H_gtr_dp","P_gtr_beta","P_gtr_xp","P_gtr_yp","P_gtr_p","P_gtr_dp"]
#    HMS_Data_Header = ["H_gtr_beta","H_gtr_xp","H_gtr_yp","H_gtr_dp"]
#    SHMS_Data_Header = ["P_gtr_beta","P_gtr_xp","P_gtr_yp","P_gtr_p","P_gtr_dp"]

    All_Data_Header = ["H_hod_goodscinhit", "H_hod_goodstarttime", "H_hod_betanotrack", "H_gtr_xp", "H_gtr_yp", "H_gtr_dp", "H_hgcer_npeSum", "H_cal_etotnorm", "H_dc_ntrack", "P_hod_goodscinhit", "P_hod_goodstarttime", "P_hod_betanotrack", "P_gtr_xp", "P_gtr_yp", "P_gtr_dp", "P_hgcer_npeSum", "P_aero_npeSum", "P_cal_etotnorm", "P_dc_Ch1_nhit", "P_dc_Ch2_nhit",  "P_dc_x_fp", "P_dc_y_fp", "P_dc_xp_fp", "P_dc_yp_fp", "XDip", "YDip", "P_dc_ntrack"]
    SHMS_Only_Data_Header = ["P_hod_goodscinhit", "P_hod_goodstarttime", "P_hod_betanotrack", "P_gtr_xp", "P_gtr_yp", "P_gtr_dp", "P_hgcer_npeSum", "P_aero_npeSum", "P_cal_etotnorm", "P_dc_Ch1_nhit", "P_dc_Ch2_nhit", "P_dc_x_fp", "P_dc_y_fp", "P_dc_xp_fp", "P_dc_yp_fp", "XDip", "YDip", "P_dc_ntrack"]
    Raw_Ctime_Data_Header = ["C_raw_ctime", "P_hod_goodscinhit", "H_hod_goodscinhit"]
    COIN_Data_Header = ["H_hod_goodscinhit", "H_hod_goodstarttime", "H_hod_betanotrack", "H_gtr_xp", "H_gtr_yp", "H_gtr_dp", "H_hgcer_npeSum", "H_cal_etotnorm", "H_dc_ntrack", "P_hod_goodscinhit", "P_hod_goodstarttime", "P_hod_betanotrack", "P_gtr_xp", "P_gtr_yp", "P_gtr_dp", "P_hgcer_npeSum", "P_aero_npeSum", "P_cal_etotnorm", "P_dc_Ch1_nhit", "P_dc_Ch2_nhit", "P_dc_x_fp", "P_dc_y_fp", "P_dc_xp_fp", "P_dc_yp_fp", "XDip", "YDip", "P_dc_ntrack"]
    HMS_SHMS_Track1_Data_Header = ["H_hod_goodscinhit", "H_hod_goodstarttime", "H_hod_betanotrack", "H_gtr_xp", "H_gtr_yp", "H_gtr_dp", "H_hgcer_npeSum", "H_cal_etotnorm", "H_dc_ntrack", "P_hod_goodscinhit", "P_hod_goodstarttime", "P_hod_betanotrack", "P_gtr_xp", "P_gtr_yp", "P_gtr_dp", "P_hgcer_npeSum", "P_aero_npeSum", "P_cal_etotnorm", "P_dc_Ch1_nhit", "P_dc_Ch2_nhit", "P_dc_x_fp", "P_dc_y_fp", "P_dc_xp_fp", "P_dc_yp_fp", "XDip", "YDip", "P_dc_ntrack"]
    HMS_SHMS_Trackn_Data_Header = ["H_hod_goodscinhit", "H_hod_goodstarttime", "H_hod_betanotrack", "H_gtr_xp", "H_gtr_yp", "H_gtr_dp", "H_hgcer_npeSum", "H_cal_etotnorm", "H_dc_ntrack", "P_hod_goodscinhit", "P_hod_goodstarttime", "P_hod_betanotrack", "P_gtr_xp", "P_gtr_yp", "P_gtr_dp", "P_hgcer_npeSum", "P_aero_npeSum", "P_cal_etotnorm", "P_dc_Ch1_nhit", "P_dc_Ch2_nhit", "P_dc_x_fp", "P_dc_y_fp", "P_dc_xp_fp", "P_dc_yp_fp", "XDip", "YDip", "P_dc_ntrack"]
    

    data = {} # Create an empty dictionary
    for d in (All_Events_Data, SHMS_Only_Events_Data, Raw_Ctime_Events_Data, COIN_Events_Data, HMS_SHMS_Track1_Events_Data, HMS_SHMS_Trackn_Events_Data): # Convert individual dictionaries into a "dict of dicts"
        data.update(d) # For every dictionary we give above, add its keys to the new dict
        data_keys = list(data.keys()) # Create a list of all the keys in all dicts added above, each is an array of data

    for i in range (0, len(data_keys)):
        # Set the headers for our data frame
        if ("All_" in data_keys[i]):
            DFHeader=list(All_Data_Header)
        elif("SHMS_Only_Events" in data_keys[i]):
            DFHeader=list(SHMS_Only_Data_Header)
        elif("Raw_Ctime_Events" in data_keys[i]):
            DFHeader=list(Raw_Ctime_Data_Header)
        elif("COIN_Events" in data_keys[i]):
            DFHeader=list(COIN_Data_Header)
        elif("Track1_" in data_keys[i]):
            DFHeader=list(HMS_SHMS_Track1_Data_Header)
        elif("Trackn_" in data_keys[i]):
            DFHeader=list(HMS_SHMS_Trackn_Data_Header)
        else:
            continue
        if (i == 0): # For the first case, start writing to file
            pd.DataFrame(data.get(data_keys[i]), columns = DFHeader, index = None).to_root("%s/%s_%s_track_Data.root" % (OUTPATH, runNum, MaxEvent), key ="%s" % data_keys[i])
        if (i != 0): # For any but the first case, append it to our file
            pd.DataFrame(data.get(data_keys[i]), columns = DFHeader, index = None).to_root("%s/%s_%s_track_Data.root" % (OUTPATH, runNum, MaxEvent), key ="%s" % data_keys[i], mode ='a') 
                    
if __name__ == '__main__':
    main() 
