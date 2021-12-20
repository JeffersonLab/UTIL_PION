#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2021-12-15 06:33:23 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

# 15/01/21 - Stephen Kay, University of Regina
# 21/06/21 - Edited By - Muhammad Junaid, University of Regina, Canada

# Python version of the pion analysis script. Now utilises uproot to select event of each type and writes them to a root file
# Intention is to apply PID/selection cutting here and plot in a separate script
# Python should allow for easier reading of databases storing timing offsets e.t.c.
# 27/04/21 - Updated to use new hcana variables, old determinations removed

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


##################################################################################################################################################

# Check the number of arguments provided to the script
if len(sys.argv)-1!=4:
    print("!!!!! ERROR !!!!!\n Expected 4 arguments\n Usage is with - ROOTfilePrefix RunNumber MaxEvents spec \n!!!!! ERROR !!!!!")
    sys.exit(1)

##################################################################################################################################################

# Input params - run number and max number of events
ROOTPrefix = sys.argv[1]
runNum = sys.argv[2]
MaxEvent = sys.argv[3]
spec = sys.argv[4]

spec = spec.upper()

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''

# Import package for cuts
import ltsep as lt 

# Add this to all files for more dynamic pathing
USER =  lt.SetPath(os.path.realpath(__file__)).getPath("USER") # Grab user info for file finding
HOST = lt.SetPath(os.path.realpath(__file__)).getPath("HOST")
REPLAYPATH = lt.SetPath(os.path.realpath(__file__)).getPath("REPLAYPATH")
UTILPATH = lt.SetPath(os.path.realpath(__file__)).getPath("UTILPATH")
ANATYPE=lt.SetPath(os.path.realpath(__file__)).getPath("ANATYPE")

################################################################################################################################################

# Add more path setting as needed in a similar manner
OUTPATH = "%s/OUTPUT/Analysis/HeeP" % UTILPATH        # Output folder location
CUTPATH = "%s/DB/CUTS" % UTILPATH

################################################################################################################################################
'''
Check that root/output paths and files exist for use
'''

# Construct the name of the rootfile based upon the info we provided
rootName = "%s/ROOTfiles/Analysis/HeeP/%s_%s_%s.root" % (UTILPATH, ROOTPrefix, runNum, MaxEvent)     # Input file location and variables taking
print ("Attempting to process %s" %(rootName))
lt.SetPath(os.path.realpath(__file__)).checkDir(OUTPATH)
lt.SetPath(os.path.realpath(__file__)).checkFile(rootName)
print("Output path checks out, outputting to %s" % (OUTPATH))

###############################################################################################################################################

# Read stuff from the main event tree
e_tree = up.open(rootName)["T"]

if spec == "HMS":
    # HMS info
    H_hod_goodscinhit = e_tree.array("H.hod.goodscinhit")            #
    H_hod_goodstarttime = e_tree.array("H.hod.goodstarttime")        #
    H_gtr_beta = e_tree.array("H.gtr.beta")                          # Beta is velocity of particle between pairs of hodoscopes
    H_gtr_xp = e_tree.array("H.gtr.th")                              # xpfp -> Theta
    H_gtr_yp = e_tree.array("H.gtr.ph")                              # ypfp -> Phi
    H_gtr_dp = e_tree.array("H.gtr.dp")                              # dp is Delta
    H_gtr_p = e_tree.array("H.gtr.p")                                # 
    H_cal_etotnorm = e_tree.array("H.cal.etotnorm")                  #
    H_cal_etottracknorm = e_tree.array("H.cal.etottracknorm")        #
    H_cer_npeSum = e_tree.array("H.cer.npeSum")                      #
    H_RF_Dist = e_tree.array("RFTime.HMS_RFtimeDist")                #                   
    H_W = e_tree.array("H.kin.primary.W")                              # JM added in, 2021/12/11

if spec == "SHMS":    
    # SHMS info
    P_hod_goodscinhit = e_tree.array("P.hod.goodscinhit")            #
    P_hod_goodstarttime = e_tree.array("P.hod.goodstarttime")        #
    P_gtr_beta = e_tree.array("P.gtr.beta")                          # Beta is velocity of particle between pairs of hodoscopes
    P_gtr_xp = e_tree.array("P.gtr.th")                              # xpfp -> Theta
    P_gtr_yp = e_tree.array("P.gtr.ph")                              # ypfp -> Phi
    P_gtr_dp = e_tree.array("P.gtr.dp")                              # dp is Delta
    P_gtr_p = e_tree.array("P.gtr.p")                                # 
    P_cal_etotnorm = e_tree.array("P.cal.etotnorm")                  #
    P_cal_etottracknorm = e_tree.array("P.cal.etottracknorm")        #
    P_aero_npeSum = e_tree.array("P.aero.npeSum")                    #
    P_aero_xAtAero = e_tree.array("P.aero.xAtAero")                  #
    P_aero_yAtAero = e_tree.array("P.aero.yAtAero")                  #
    P_hgcer_npeSum = e_tree.array("P.hgcer.npeSum")                  #
    P_hgcer_xAtCer = e_tree.array("P.hgcer.xAtCer")                  #
    P_hgcer_yAtCer = e_tree.array("P.hgcer.yAtCer")                  #
    P_ngcer_npeSum = e_tree.array("P.ngcer.npeSum")                  #
    P_ngcer_xAtCer = e_tree.array("P.ngcer.xAtCer")                  #
    P_ngcer_yAtCer = e_tree.array("P.ngcer.yAtCer")                  #
    P_RF_Dist = e_tree.array("RFTime.SHMS_RFtimeDist")               #
    emiss = e_tree.array("P.kin.secondary.emiss")                   
    pmiss = e_tree.array("P.kin.secondary.pmiss")                   
    MMpi = e_tree.array("P.kin.secondary.MMpi")                      
    W = e_tree.array("P.kin.primary.W")                              
    pmiss_x = e_tree.array("P.kin.secondary.pmiss_x")                
    pmiss_y = e_tree.array("P.kin.secondary.pmiss_y")                
    pmiss_z = e_tree.array("P.kin.secondary.pmiss_z")                


##############################################################################################################################################

# Defining path for cut file  
if spec == "HMS":
    fout = '%s/DB/CUTS/run_type/hSing_prod.cuts' % UTILPATH
if spec == "SHMS":
    fout = '%s/DB/CUTS/run_type/pSing_prod.cuts' % UTILPATH

#################################################################################################################################################################

# defining Cuts
if spec == "HMS":
    cuts = ["sing_ee_cut_all_noRF"]
if spec == "SHMS":
    cuts = ["sing_ee_cut_ngcer_all_noRF"]

def make_cutDict(cuts,fout,runNum,CURRENT_ENV,DEBUG=False):
    '''
    This method calls several methods in kaonlt package. It is required to create properly formated
    dictionaries. The evaluation must be in the analysis script because the analysis variables (i.e. the
    leaves of interest) are not defined in the kaonlt package. This makes the system more flexible
    overall, but a bit more cumbersome in the analysis script. Perhaps one day a better solution will be
    implimented.
    '''

    # read in cuts file and make dictionary
    importDict = lt.SetCuts(CURRENT_ENV).importDict(cuts,fout,runNum,False)
    for i,cut in enumerate(cuts):
        x = lt.SetCuts(CURRENT_ENV,importDict).booleanDict(cut)
        print("\n%s" % cut)
        print(x, "\n")
        if i == 0:
            inputDict = {}
        cutDict = lt.SetCuts(CURRENT_ENV,importDict).readDict(cut,inputDict)
        for j,val in enumerate(x):
            cutDict = lt.SetCuts(CURRENT_ENV,importDict).evalDict(cut,eval(x[j]),cutDict)
    return lt.SetCuts(CURRENT_ENV,cutDict)

c = make_cutDict(cuts,fout,runNum,os.path.realpath(__file__))

#################################################################################################################################################################

def sing():

    if spec == "HMS":
        # Define the array of arrays containing the relevant HMS and SHMS info                              
        
        NoCut_SING = [H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_gtr_p, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, H_RF_Dist]
        
        Uncut_SING = [(H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_gtr_p, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, H_RF_Dist) for (H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_gtr_p, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, H_RF_Dist) in zip(*NoCut_SING)]
        
        # Create array of arrays of pions after cuts, all events
        
        Cut_SING_tmp = NoCut_SING
        Cut_SING_all_tmp = []
        
        for arr in Cut_SING_tmp:
            Cut_SING_all_tmp.append(c.add_cut(arr, "sing_ee_cut_all_noRF"))
            
        Cut_SING_all = [(H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_gtr_p, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, H_RF_Dist) for (H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_gtr_p, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, H_RF_Dist) in zip(*Cut_SING_all_tmp)]
        
        SING = {
            "Uncut_Events" : Uncut_SING,
            "Cut_Events_All" : Cut_SING_all,
        }

    if spec == "SHMS":
        # Define the array of arrays containing the relevant HMS and SHMS info                              
        
        NoCut_SING = [pmiss_z, pmiss_y, pmiss_x, W, MMpi, pmiss, emiss, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, P_ngcer_npeSum, P_ngcer_xAtCer, P_ngcer_yAtCer, P_RF_Dist]
        
        Uncut_SING = [(pmiss_z, pmiss_y, pmiss_x, W, MMpi, pmiss, emiss, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, P_ngcer_npeSum, P_ngcer_xAtCer, P_ngcer_yAtCer, P_RF_Dist) for (pmiss_z, pmiss_y, pmiss_x, W, MMpi, pmiss, emiss, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, P_ngcer_npeSum, P_ngcer_xAtCer, P_ngcer_yAtCer, P_RF_Dist) in zip(*NoCut_SING)]
        
        # Create array of arrays of pions after cuts, all events
        
        Cut_SING_tmp = NoCut_SING
        Cut_SING_all_tmp = []
        
        for arr in Cut_SING_tmp:
            Cut_SING_all_tmp.append(c.add_cut(arr, "sing_ee_cut_ngcer_all_noRF"))
            
        Cut_SING_all = [(pmiss_z, pmiss_y, pmiss_x, W, MMpi, pmiss, emiss, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, P_ngcer_npeSum, P_ngcer_xAtCer, P_ngcer_yAtCer, P_RF_Dist) for (pmiss_z, pmiss_y, pmiss_x, W, MMpi, pmiss, emiss, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, P_ngcer_npeSum, P_ngcer_xAtCer, P_ngcer_yAtCer, P_RF_Dist) in zip(*Cut_SING_all_tmp)]
        
        SING = {
            "Uncut_Events" : Uncut_SING,
            "Cut_Events_All" : Cut_SING_all,
        }
                
    return SING
        
##################################################################################################################################################################

def main():
    SING_Data = sing()

    # This is just the list of branches we use from the initial root file for each dict
    # I don't like re-defining this here as it's very prone to errors if you included (or removed something) earlier but didn't modify it here
    # Should base the branches to include based on some list and just repeat the list here (or call it again directly below)

    if spec == "HMS":
        SING_Data_Header = ["H_gtr_beta","H_gtr_xp","H_gtr_yp","H_gtr_dp", "H_gtr_p","H_hod_goodscinhit","H_hod_goodstarttime","H_cal_etotnorm","H_cal_etottracknorm","H_cer_npeSum","H_RF_Dist"]
    if spec == "SHMS":
        SING_Data_Header = ["pmiss_z","pmiss_y","pmiss_x", "W","MMpi","pmiss","emiss","P_gtr_beta","P_gtr_xp","P_gtr_yp","P_gtr_p","P_gtr_dp","P_hod_goodscinhit","P_hod_goodstarttime","P_cal_etotnorm","P_cal_etottracknorm","P_aero_npeSum","P_aero_xAtAero","P_aero_yAtAero","P_hgcer_npeSum","P_hgcer_xAtCer","P_hgcer_yAtCer", "P_ngcer_npeSum", "P_ngcer_xAtCer", "P_ngcer_yAtCer","P_RF_Dist"]
        
    # Need to create a dict for all the branches we grab                                                
    data = {}
    data.update(SING_Data)
    data_keys = list(data.keys()) # Create a list of all the keys in all dicts added above, each is an array of data                                                                                       

    for i in range (0, len(data_keys)):
        if("Events" in data_keys[i]):
            DFHeader=list(SING_Data_Header)
        else:
            continue
            # Uncomment the line below if you want .csv file output, WARNING the files can be very large and take a long time to process!                                                                      
            #pd.DataFrame(data.get(data_keys[i])).to_csv("%s/%s_%s.csv" % (OUTPATH, data_keys[i], runNum), header=DFHeader, index=False) # Convert array to panda dataframe and write to csv with correct header                                                                                                      
        if (i == 0):
            pd.DataFrame(data.get(data_keys[i]), columns = DFHeader, index = None).to_root("%s/%s_%s_%s_An\
alysed_Data.root" % (OUTPATH, spec, runNum, MaxEvent), key ="%s" % data_keys[i])
        elif (i != 0):
            pd.DataFrame(data.get(data_keys[i]), columns = DFHeader, index = None).to_root("%s/%s_%s_%s_An\
alysed_Data.root" % (OUTPATH, spec, runNum, MaxEvent), key ="%s" % data_keys[i], mode ='a')

if __name__ == '__main__':
    main()
print ("Processing Complete")
