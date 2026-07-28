#! /usr/bin/python
#
# Description: Grabs lumi data from corresponding csv depending on run setting. Then plots the yields and creates a comprehensive table.
# Variables calculated: current, rate_HMS, rate_SHMS, sent_edtm_PS, uncern_HMS_evts_scaler, uncern_SHMS_evts_scaler, uncern_HMS_evts_notrack, uncern_SHMS_evts_notrack, uncern_HMS_evts_track, uncern_SHMS_evts_track
# ================================================================
# Time-stamp: "2023-08-31 16:43:47 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
# Modified to handle coin lumi case 
# by heinricn - 2024/09/19
#
import numpy as np  # type: ignore
import pandas as pd # type: ignore
import sys, os, subprocess, math # import this line so you can use os.environ.get

import matplotlib as mpl                # type: ignore # define this line before .pyplot 
if os.environ.get('DISPLAY', '')=='':   # Added to "fix" the "no display name and no $DISPLAY environment variable" error
    print('no display found. Using non-interactive Agg backend') #
    mpl.use('Agg')                      # causes output in terminal "No handles with labels found to put in legend"

import matplotlib.pyplot as plt # type: ignore
from scipy.optimize import curve_fit # type: ignore
#import statsmodels.api as sm # type: ignore
from csv import DictReader



################################################################################################################################################
'''
ltsep package import and pathing definitions
'''

# Import package for cuts
from ltsep import Root # type: ignore

lt=Root(os.path.realpath(__file__))

# Add this to all files for more dynamic pathing
USER=lt.USER # Grab user info for file finding
HOST=lt.HOST
REPLAYPATH=lt.REPLAYPATH
UTILPATH=lt.UTILPATH
SCRIPTPATH=lt.SCRIPTPATH
ANATYPE=lt.ANATYPE

################################################################################################################################################
'''
Grab proper lumi data file based of input given then define input and output file names as well as the target
'''

inp_name = sys.argv[1]

sys.path.insert(0,"%s/luminosity/src/%sLT" % (SCRIPTPATH,ANATYPE))
import data_path # type: ignore

data_path = data_path.get_file(inp_name,SCRIPTPATH)
target = data_path[0]
inp_f = data_path[1]
out_f = data_path[2]

#assume both HMS and SHMS data in file
NO_HMS = False
NO_SHMS = False

################################################################################################################################################

print("\nRunning as %s on %s, hallc_replay_lt path assumed as %s" % (USER, HOST, REPLAYPATH))

# Converts csv data to dataframe
try:
    lumi_data = pd.read_csv(inp_f)
    print(inp_f)
    print(lumi_data.keys())
    print(lumi_data["tot_events"])
except IOError:
    print("Error: %s does not appear to exist." % inp_f)
    sys.exit(0)

################################################################################################################################################

# Remove runs that have zero beam on time to avoid error
lumi_data = lumi_data[lumi_data["time"] > 0.0].reset_index(drop=True)
# Removing runs that are less than a minute
# Why? - NH
#lumi_data = lumi_data[ (lumi_data['time'] >= 60.0) ].reset_index(drop=True)

def removeRun(runNum):
    '''
    Removes runs from DF and subsequently will not be plotted or included in yield csv output
    '''
    global lumi_data
    if (lumi_data['run number'] == runNum).any():
        print("Before",lumi_data["run number"].values)
        print("Removing run {} from lumi_data...".format(runNum))
        lumi_data = lumi_data[lumi_data['run number'] != runNum].reset_index(drop=True)
        print("After",lumi_data["run number"].values)


################################################################################################################################################

# Convert to dict for proper formatting when eventually merging dictionaries
lumi_data = dict(lumi_data)
    
################################################################################################################################################

# Define number of runs to analyze
numRuns = len(lumi_data["run number"])

def makeList(lumi_input):
    '''
    Takes input list and converts to numpy (with NaN converted to zeros) so it can be used mathematically
    '''
    new_lst = [lumi_data[lumi_input][i] for i,evts in enumerate(lumi_data["run number"])]
    new_lst = np.asarray(pd.Series(new_lst).fillna(0)) # changes NaN to zeros and convert to numpy
    return new_lst

################################################################################################################################################

def calc_yield():
    '''
    Creates a new dictionary with yield calculations. The relative yield is defined relative to the maximum current.
    '''

    # Create dictionary for calculations that were not calculated in previous scripts.
    yield_dict = {
        "current" : makeList("charge")/makeList("time"),
        
        "rate_HMS" : makeList("HMSTRIG_scaler")/makeList("time"),
        "rate_SHMS" : makeList("SHMSTRIG_scaler")/makeList("time"),
        "rate_COIN" : makeList("COINTRIG_scaler")/makeList("time"),
        "rate_TOTAL" : (makeList("COINTRIG_scaler") + makeList("SHMSTRIG_scaler") + makeList("HMSTRIG_scaler"))/makeList("time"), #for ELT live time
        
        "rate_HMS3-4" : makeList("HMS_3-4Trig_scaler")/makeList("time"),
        "rate_SHMS3-4" : makeList("SHMS_3-4Trig_scaler")/makeList("time"),
        "uncern_HMS3-4_evts_scaler" : np.sqrt(makeList("HMS_3-4Trig_scaler"))/makeList("HMS_3-4Trig_scaler"),
        "uncern_SHMS3-4_evts_scaler" : np.sqrt(makeList("SHMS_3-4Trig_scaler"))/makeList("SHMS_3-4Trig_scaler"),
        
        "SHMS3of4ELT" : makeList("SHMS3of4ELT"),
        "HMS3of4ELT" : makeList("HMS3of4ELT"),
        "SHMS3of4ELT_err" : makeList("SHMS3of4ELT_err"),
        "HMS3of4ELT_err" : makeList("HMS3of4ELT_err"),
        
        "CoinBlocking_Corr": makeList("CoinBlocking_Corr"),
        "Uncer_CoinBlocking_Corr": np.sqrt(makeList("CoinBlocking_Corr")*(1- makeList("CoinBlocking_Corr"))/makeList("CoinBlocking_Total")),
        "SysUncer_CoinBlocking_Corr_up":makeList("CoinBlocking_SysOut")/makeList("CoinBlocking_Total") - makeList("CoinBlocking_Corr"),
        "SysUncer_CoinBlocking_Corr_down":makeList("CoinBlocking_SysIn")/makeList("CoinBlocking_Total") - makeList("CoinBlocking_Corr"),
        "SysUncer_CoinBlocking_Corr":(makeList("CoinBlocking_SysOut")/makeList("CoinBlocking_Total") - makeList("CoinBlocking_SysIn")/makeList("CoinBlocking_Total"))/2,

        #"CPULT_phys" : (makeList("HMSTRIG_cut")*makeList("HMS_PS")+makeList("SHMSTRIG_cut")*makeList("SHMS_PS"))*makeList("CPULT_scaler"),
        "CPULT_phys" : makeList("CPULT_scaler"),

        "uncern_CPULT_phys" : makeList("CPULT_scaler_uncern"),
        
        #"ELT_3of4_SHMS" : makeList("SHMS3of4ELT"),
        
        "uncern_HMS_evts_scaler" : np.sqrt(makeList("HMSTRIG_scaler"))/makeList("HMSTRIG_scaler"),
        
        "uncern_SHMS_evts_scaler" : np.sqrt(makeList("SHMSTRIG_scaler"))/makeList("SHMSTRIG_scaler"),

        "uncern_COIN_evts_scaler" : np.sqrt(makeList("COINTRIG_scaler"))/makeList("COINTRIG_scaler"),

        "uncern_HMS_evts_notrack" : np.sqrt(abs(makeList("h_int_etotnorm_evts")))/makeList("h_int_etotnorm_evts"),

        "uncern_SHMS_evts_notrack" : np.sqrt(abs(makeList("p_int_etotnorm_evts")))/makeList("p_int_etotnorm_evts"),

        "uncern_HMS_evts_track" : np.sqrt(abs(makeList("h_int_etottracknorm_evts")))/makeList("h_int_etottracknorm_evts"),

        "uncern_SHMS_evts_track" : np.sqrt(abs(makeList("p_int_etottracknorm_evts")))/makeList("p_int_etottracknorm_evts"),
        
        "uncern_COIN_evts_notrack" : np.sqrt(abs((makeList("c_int_noTrack_events")/(makeList("c_int_noTrack_events")**2) + makeList("c_int_noTrack_Rand_events")/(makeList("c_int_noTrack_Rand_events")**2))/(makeList("c_int_noTrack_events") - makeList("c_int_noTrack_Rand_events")))),
        
        "uncern_COIN_evts_track" : np.sqrt(abs((makeList("c_int_Track_events")/(makeList("c_int_Track_events")**2) + makeList("c_int_Track_Rand_events")/(makeList("c_int_Track_Rand_events")**2))/(makeList("c_int_Track_events") - makeList("c_int_Track_Rand_events"))))
    }            
    
    # defult case, setting to 2021 for now - NH
    slope = 5513.79
    uncern_slope = 7.6
    intercept = 250429
    uncern_intercept = 271
    #2022 data or 2021 data
    if "2022" in inp_name:
        slope = 5542.0
        uncern_slope = 11.5
        intercept = 250029
        uncern_intercept = 270
    else: 
        if "2021" in inp_name:
            slope = 5513.79
            uncern_slope = 7.6
            intercept = 250429
            uncern_intercept = 271
        else: pass
            
    uncern_charge = np.sqrt(((((makeList("charge")**2)*((uncern_slope**2))+ (makeList("time")**2)+uncern_intercept**2))/slope**2))
    yield_dict.update({"uncern_charge" : uncern_charge})
        
    rateHMSSlope = 4.5/(10**7) #adhoc rate correction for 0.04%/kHz slope
    rateHMSSlopeUncer = 0.22/(10**7)  #statistical uncertainty
    
    rateSHMSSlope = 4.5/(10**7) #adhoc rate correction for 0.04%/kHz slope
    rateSHMSSlopeUncer = 0.22/(10**7)  #statistical uncertainty
    
    rateCOINslope = 0.030298/(10**3) #based off off pass 22 coin rate slope. 3.1%/kHz sl
    rateCOINslopeUncer = (0.000535 + 0.00402)/(10**3)

    #rateHMSCorr = 1 - (yield_dict["rate_HMS3-4"]*(rateHMSSlope))         
    #uncer_rateHMSCorr = (yield_dict["rate_HMS3-4"]*rateHMSSlopeUncer)**2 + (rateHMSSlope*yield_dict["uncern_HMS3-4_evts_scaler"])**2
    #rateSHMSCorr = 1 - (yield_dict["rate_SHMS3-4"]*(rateSHMSSlope))         
    #uncer_rateSHMSCorr = (yield_dict["rate_SHMS3-4"]*rateSHMSSlopeUncer)**2 + (rateSHMSSlope*yield_dict["uncern_SHMS3-4_evts_scaler"])**2
    
    rateHMSCorr = (yield_dict["HMS3of4ELT"])         
    uncer_rateHMSCorr = yield_dict["HMS3of4ELT_err"]
    rateSHMSCorr = (yield_dict["SHMS3of4ELT"])         
    uncer_rateSHMSCorr = yield_dict["SHMS3of4ELT_err"]
    #rateCOINCorr = 1-(rateCOINslope*yield_dict["rate_COIN"])
    #rateCOINCorrUncer = (rateCOINslopeUncer*yield_dict["rate_COIN"])
    rateCOINCorr = 1
    rateCOINCorrUncer = 0 

    #temporary to test with just boiling correction
    #rateHMSCorr = rateHMSCorr/rateHMSCorr
    #rateSHMSCorr = rateSHMSCorr/rateSHMSCorr
    #uncer_rateHMSCorr = rateHMSCorr*0
    #uncer_rateSHMSCorr = rateSHMSCorr*0
    
    yield_dict.update({"rateHMSCorr": rateHMSCorr})
    yield_dict.update({"uncer_rateHMSCorr": uncer_rateHMSCorr})
    yield_dict.update({"rateSHMSCorr": rateSHMSCorr})
    yield_dict.update({"uncer_rateSHMSCorr": uncer_rateSHMSCorr})
    yield_dict.update({"rateCOINCorr": rateCOINCorr})
    yield_dict.update({"rateCOINCorrUncer": rateCOINCorrUncer})
        
    boilSlope = 2.8/10000 #boiling correction for 2.8%/100uA slope
    boilSlopeUncerStat = 0.17/10000 #statistical + systematic uncertainty
    boilSlopeUncerSys = 0.32/10000 #statistical + systematic uncertainty
    
    boilingCorr = 1 - (yield_dict["current"]*makeList("curr_corr")*(boilSlope)) 
    uncer_Current =  (((uncern_slope**2)+ (uncern_intercept**2/yield_dict["current"]))/slope**2)
    uncer_boilingCorrStat =  (yield_dict["current"]*makeList("curr_corr")*boilSlopeUncerStat)**2 + (boilSlope*uncer_Current)**2
    uncer_boilingCorrSys =  (yield_dict["current"]*makeList("curr_corr")*boilSlopeUncerSys)**2 + (boilSlope*uncer_Current)**2
    uncer_boilingCorr =  uncer_boilingCorrSys + uncer_boilingCorrStat

    yield_dict.update({"boilingCorr": boilingCorr})
    yield_dict.update({"uncer_boilingCorrStat": uncer_boilingCorrStat})
    yield_dict.update({"uncer_boilingCorrSys": uncer_boilingCorrSys})
    yield_dict.update({"uncer_boilingCorr": uncer_boilingCorr})

    gwidth = 55*(10**(-9)) #gate width == 55 nanoseconds
    uncern_gwidth = 5*(10**(-9)) # 5 ns uncerntainty
        
    # Total livetime calculation
    TLT_EDTM = makeList("accp_edtm")/makeList("sent_edtm_PS") #TLT with edtm 
    
    if makeList("HMS_PS").any() != 0:
        yield_dict.update({"rate_HMS": makeList("HMSTRIG_scaler")/makeList("time")})
    else:  yield_dict.update({"rate_HMS": 0})
        
    if makeList("SHMS_PS").any() != 0:
        yield_dict.update({"rate_SHMS": makeList("SHMSTRIG_scaler")/makeList("time")})  
    else: yield_dict.update({"rate_SHMS": 0})
    
    EDT2 = gwidth*(yield_dict["rate_HMS"] + yield_dict["rate_SHMS"]) # estimation of dead time using gate width and trigger rate
    EDT2_HMS = gwidth*(yield_dict["rate_HMS"]) # estimation of dead time using gate width and trigger rate
    EDT2_SHMS = gwidth*(yield_dict["rate_SHMS"]) # estimation of dead time using gate width and trigger rate

    ELT2_HMS  = 1 - EDT2_HMS
    ELT2_SHMS = 1 - EDT2_SHMS
    uncern_ELT2_HMS = (ELT2_HMS*(1-ELT2_HMS))/(50*ELT2_HMS**2)
        
    
    TLT_ELT_SHMS = yield_dict["rateSHMSCorr"]*yield_dict["CPULT_phys"]
    TLT_ELT_HMS = yield_dict["rateHMSCorr"]*yield_dict["CPULT_phys"]
    TLT_ELT = yield_dict["rateSHMSCorr"]*yield_dict["rateHMSCorr"]*yield_dict["CPULT_phys"]*yield_dict["rateCOINCorr"]

    CPULT_TLT_SHMS = TLT_EDTM/TLT_ELT_SHMS
    CPULT_TLT_HMS = TLT_EDTM/TLT_ELT_HMS

    yield_dict.update({"CPULT_TLT_SHMS" : CPULT_TLT_SHMS})
    yield_dict.update({"CPULT_TLT_HMS" : CPULT_TLT_HMS})
    yield_dict.update({"TLT" : TLT_EDTM})
    yield_dict.update({"TLT_ELT_SHMS" : TLT_ELT_SHMS})
    yield_dict.update({"TLT_ELT_HMS" : TLT_ELT_HMS})
    yield_dict.update({"TLT_ELT" : TLT_ELT})
    yield_dict.update({"ELT2_HMS" : ELT2_HMS})
    yield_dict.update({"ELT2_SHMS" : ELT2_SHMS})
    yield_dict.update({"uncern_ELT2_HMS" : uncern_ELT2_HMS})

    #yield_dict.update({"TLT" : TLT})
    
    # Stat + Systematic errors
    uncern_TLT = np.sqrt(yield_dict["TLT"]*(1-yield_dict["TLT"])/makeList("accp_edtm"))
   
    uncern_TLT_ELT_SHMS = ((makeList("CPULT_scaler_uncern")/makeList("CPULT_scaler"))**2 + (yield_dict["uncer_rateSHMSCorr"]/yield_dict["rateSHMSCorr"])**2)
    uncern_TLT_ELT_HMS = ((makeList("CPULT_scaler_uncern")/makeList("CPULT_scaler"))**2 + (yield_dict["uncer_rateHMSCorr"]/yield_dict["rateHMSCorr"])**2)
    uncern_TLT_ELT = (makeList("CPULT_scaler_uncern")/makeList("CPULT_scaler"))**2 + (yield_dict["uncer_rateSHMSCorr"]/yield_dict["rateSHMSCorr"])**2 + (yield_dict["uncer_rateHMSCorr"]/yield_dict["rateHMSCorr"])**2

    #uncern_TLT_ELT = np.sqrt((ELT2*(1-ELT2))/100)
    
    #uncern_TLT = np.sqrt(makeList("sent_edtm_PS")*.95*.05)
    yield_dict.update({"uncern_TLT" : uncern_TLT})
    yield_dict.update({"uncern_TLT_ELT" : uncern_TLT_ELT})
    yield_dict.update({"uncern_TLT_ELT_SHMS" : uncern_TLT_ELT_SHMS })
    yield_dict.update({"uncern_TLT_ELT_HMS" : uncern_TLT_ELT_HMS })
    
        
    yield_dict.update({"rate_SHMS" : makeList("SHMSTRIG_scaler")/makeList("time")})
    
    if makeList("SHMS_PS").any() != 0:
        yield_dict.update({"sent_edtm_SHMS" : makeList("sent_edtm")/makeList("SHMS_PS")}) # get rid of divide by zero error for single arm runs
    else: yield_dict.update({"sent_edtm_SHMS" : 1})

    yield_dict.update({"uncern_SHMS_evts_scaler" : np.sqrt(makeList("SHMSTRIG_scaler"))/makeList("SHMSTRIG_scaler")})
    yield_dict.update({"uncern_SHMS_evts_notrack" : np.sqrt(abs(makeList("p_int_etotnorm_evts")))/makeList("p_int_etotnorm_evts")})
    yield_dict.update({"uncern_SHMS_evts_track" : np.sqrt(abs(makeList("p_int_etottracknorm_evts")))/makeList("p_int_etottracknorm_evts")})

    pTLT = makeList("paccp_edtm")/yield_dict["sent_edtm_SHMS"]
    yield_dict.update({"pTLT" : pTLT})        

    COIN_scaler_accp = makeList("COINTRIG_scaler")-makeList("sent_edtm")
    yield_COIN_scaler = COIN_scaler_accp/(makeList("charge")*makeList("curr_corr"))
    #yield_COIN_notrack = (makeList("c_int_noTrack_events")-(makeList("c_int_noTrack_Rand_events")/6))/(makeList("charge")*makeList("curr_corr")*yield_dict["TLT"]*yield_dict["boilingCorr"]*yield_dict["rateHMSCorr"]*yield_dict["rateSHMSCorr"])
    #yield_COIN_track = (makeList("c_int_Track_events")-(makeList("c_int_Track_Rand_events")/6))/(makeList("charge")*makeList("curr_corr")*yield_dict["TLT"]*makeList("SHMS_track")*makeList("HMS_track")*yield_dict["boilingCorr"]*yield_dict["rateHMSCorr"]*yield_dict["rateSHMSCorr"])
    
    yield_COIN_notrack = (makeList("c_int_noTrack_events")-(makeList("c_int_noTrack_Rand_events")/6))/(makeList("charge")*makeList("curr_corr")*yield_dict["TLT_ELT"]*yield_dict["boilingCorr"])
    #yield_COIN_track = (makeList("c_int_Track_events")-(makeList("c_int_Track_Rand_events")/6))/(makeList("charge")*makeList("curr_corr")*yield_dict["TLT_ELT"]*yield_dict["CoinBlocking_Corr"]*makeList("SHMS_track")*makeList("HMS_track")*yield_dict["boilingCorr"])
    yield_COIN_track = (makeList("c_int_Track_events")-(makeList("c_int_Track_Rand_events")/6))/(makeList("charge")*makeList("curr_corr")*yield_dict["TLT"]*yield_dict["CoinBlocking_Corr"]*makeList("SHMS_track")*makeList("HMS_track")*yield_dict["boilingCorr"])
    
    print("\n\n******************")
    print("yield_COIN_notrack:")
    print(yield_COIN_notrack)
    print("yield_COIN_track")
    print(yield_COIN_track)
    print("\n\n*******************")
    
    yield_dict.update({"yield_COIN_scaler" : yield_COIN_scaler})
    yield_dict.update({"yield_COIN_notrack" : yield_COIN_notrack})
    yield_dict.update({"yield_COIN_track" : yield_COIN_track})
    
    uncern_yieldRel_COIN_scaler = np.sqrt((yield_dict["uncern_COIN_evts_scaler"]/makeList("COINTRIG_scaler"))**2 + (yield_dict["uncern_charge"]/makeList("charge"))**2)
    uncern_yieldRel_COIN_notrack = np.sqrt((yield_dict["uncern_COIN_evts_notrack"]/makeList("c_int_noTrack_events"))**2 + (yield_dict["uncern_charge"]/makeList("charge"))**2 + (yield_dict["uncern_TLT_ELT"]) + yield_dict["uncer_boilingCorr"] + yield_dict["uncer_rateHMSCorr"] + yield_dict["uncer_rateSHMSCorr"])
    uncern_yieldRel_COIN_track = np.sqrt((yield_dict["uncern_COIN_evts_track"]/makeList("c_int_Track_events"))**2 + (yield_dict["uncern_charge"]/makeList("charge"))**2 + (yield_dict["uncern_TLT"]/yield_dict["TLT"])**2 + ((yield_dict["Uncer_CoinBlocking_Corr"] + yield_dict["SysUncer_CoinBlocking_Corr"])/yield_dict["CoinBlocking_Corr"])**2 + (makeList("SHMS_track_uncern")/makeList("SHMS_track"))**2 + (makeList("HMS_track_uncern")/makeList("HMS_track"))**2 + yield_dict["uncer_boilingCorr"]) 
    
    print("\n\n******************")
    print("uncern_yieldRel_COIN_track:")
    print(uncern_yieldRel_COIN_track)
    print("\n\n*******************")
    
    yield_dict.update({"uncern_yieldRel_COIN_scaler" : uncern_yieldRel_COIN_scaler})
    yield_dict.update({"uncern_yieldRel_COIN_notrack" : uncern_yieldRel_COIN_notrack})
    yield_dict.update({"uncern_yieldRel_COIN_track" : uncern_yieldRel_COIN_track})
    
    # Accepted scalers; removed EDTM subtraction. There are EDTM in here but it won't be this value.
    #SHMS_scaler_accp = makeList("SHMSTRIG_scaler") #-makeList("sent_edtm_PS")
    #SHMS_scaler_accp = makeList("SHMSTRIG_scaler")-yield_dict["sent_edtm_SHMS"]
    SHMS_scaler_accp = makeList("SHMSTRIG_scaler")-makeList("sent_edtm") # Returns reasonable yields
    yield_dict.update({"SHMS_scaler_accp" : SHMS_scaler_accp})
    if(makeList("SHMS_PS")[0] == None or makeList("SHMS_PS")[0] == -1 or makeList("SHMS_PS").any() == 0):
        NO_SHMS = True
        yield_SHMS_scaler = []
        yield_SHMS_notrack = []
        yield_SHMS_track = []
        yield_SHMS_CPULT_notrack = []
        yield_SHMS_CPULT_track = []
        for i,curr in enumerate(yield_dict["current"]):
            yield_SHMS_scaler.append(1)
            yield_SHMS_notrack.append(1)
            yield_SHMS_track.append(1)
            yield_SHMS_CPULT_notrack.append(1)
            yield_SHMS_CPULT_track.append(1)
    else:             
            yield_SHMS_scaler = (yield_dict["SHMS_scaler_accp"])/(makeList("charge")*makeList("curr_corr")*yield_dict["ELT2_SHMS"]*yield_dict["boilingCorr"]*yield_dict["rateSHMSCorr"])
            yield_SHMS_notrack = (makeList("p_int_etotnorm_evts")*makeList("SHMS_PS"))/(makeList("charge")*makeList("curr_corr")*yield_dict["TLT_ELT"]*yield_dict["boilingCorr"]*yield_dict["rateSHMSCorr"])
            yield_SHMS_track = (makeList("p_int_etottracknorm_evts")*makeList("SHMS_PS"))/(makeList("charge")*makeList("curr_corr")*yield_dict["TLT_ELT"]*makeList("SHMS_track")*yield_dict["boilingCorr"]*yield_dict["rateSHMSCorr"])
            yield_SHMS_CPULT_notrack = (makeList("p_int_etotnorm_evts")*makeList("SHMS_PS"))/(makeList("charge")*makeList("curr_corr")*yield_dict["CPULT_phys"])
            yield_SHMS_CPULT_track = (makeList("p_int_etottracknorm_evts")*makeList("SHMS_PS"))/(makeList("charge")*makeList("curr_corr")*yield_dict["CPULT_phys"]*makeList("SHMS_track"))
            
            
    yield_dict.update({"yield_SHMS_scaler" : yield_SHMS_scaler})
    yield_dict.update({"yield_SHMS_notrack" : yield_SHMS_notrack})
    yield_dict.update({"yield_SHMS_track" : yield_SHMS_track})
    yield_dict.update({"yield_SHMS_CPULT_notrack" : yield_SHMS_CPULT_notrack})
    yield_dict.update({"yield_SHMS_CPULT_track" : yield_SHMS_CPULT_track})

    #for i,curr in enumerate(yield_dict["current"]):
    #    if abs((TLT_EDTM)[i] - (TLT_ELT)[i]) > ((abs(uncern_TLT_ELT)[i]) + abs((uncern_TLT)[i])) and TLT_ELT[i] > TLT_EDTM[i]:
    #        
    #        yield_dict["yield_SHMS_notrack"][i] = (makeList("p_int_etotnorm_evts")[i]*makeList("SHMS_PS")[i])/(makeList("charge")[i]*makeList("curr_corr")[i]*yield_dict["TLT_ELT"][i])
    #        yield_dict["yield_SHMS_track"][i] = (makeList("p_int_etottracknorm_evts")[i]*makeList("SHMS_PS")[i])/(makeList("charge")[i]*makeList("curr_corr")[i]*yield_dict["TLT_ELT"][i]*makeList("SHMS_track")[i])
    #    else: pass  

    yield_dict.update({"rate_HMS" : makeList("HMSTRIG_scaler")/makeList("time")})
    
    if makeList("HMS_PS").any() != 0:
        yield_dict.update({"sent_edtm_HMS" : makeList("sent_edtm")/makeList("HMS_PS")})
    else: yield_dict.update({"sent_edtm_HMS" : 1})
    
    
    yield_dict.update({"uncern_HMS_evts_scaler" : np.sqrt(makeList("HMSTRIG_scaler"))/makeList("HMSTRIG_scaler")})
    yield_dict.update({"uncern_HMS_evts_notrack" : np.sqrt(abs(makeList("h_int_etotnorm_evts")))/makeList("h_int_etotnorm_evts")})
    yield_dict.update({"uncern_HMS_evts_track" : np.sqrt(abs(makeList("h_int_etottracknorm_evts")))/makeList("h_int_etottracknorm_evts")})

    hTLT = makeList("haccp_edtm")/yield_dict["sent_edtm_HMS"]
    yield_dict.update({"hTLT" : hTLT})

    # Accepted scalers; removed EDTM subtraction. There are EDTM in here but it won't be this value.
    #HMS_scaler_accp = makeList("HMSTRIG_scaler") #-makeList("sent_edtm_PS")
    #HMS_scaler_accp = makeList("HMSTRIG_scaler")-yield_dict["sent_edtm_HMS"]
    HMS_scaler_accp = makeList("HMSTRIG_scaler")-makeList("sent_edtm") # Returns reasonable yields    
    yield_dict.update({"HMS_scaler_accp" : HMS_scaler_accp})

    # Calculate yield values
    if(makeList("HMS_PS")[0] == None or makeList("HMS_PS")[0] == -1 or makeList("HMS_PS").any() == 0):
        NO_HMS = True
        yield_HMS_scaler = []
        yield_HMS_notrack = []
        yield_HMS_track = []
        yield_HMS_CPULT_notrack = []
        yield_HMS_CPULT_track = []
        for i,curr in enumerate(yield_dict["current"]):
            yield_HMS_scaler.append(1)
            yield_HMS_notrack.append(1)
            yield_HMS_track.append(1)
            yield_HMS_CPULT_notrack.append(1)
            yield_HMS_CPULT_track.append(1)
    else:

            yield_HMS_scaler = (yield_dict["HMS_scaler_accp"])/(makeList("charge")*makeList("curr_corr")*yield_dict["ELT2_HMS"]*yield_dict["boilingCorr"]*yield_dict["rateHMSCorr"])
            yield_HMS_notrack = (makeList("h_int_etotnorm_evts")*makeList("HMS_PS"))/(makeList("charge")*makeList("curr_corr")*yield_dict["TLT_ELT"]*yield_dict["boilingCorr"]*yield_dict["rateHMSCorr"])
            yield_HMS_track = (makeList("h_int_etottracknorm_evts")*makeList("HMS_PS"))/(makeList("charge")*makeList("curr_corr")*yield_dict["TLT_ELT"]*makeList("HMS_track")*yield_dict["boilingCorr"]*yield_dict["rateHMSCorr"])
            yield_HMS_CPULT_notrack = (makeList("h_int_etotnorm_evts")*makeList("HMS_PS"))/(makeList("charge")*makeList("curr_corr")*yield_dict["CPULT_phys"])
            yield_HMS_CPULT_track = (makeList("h_int_etottracknorm_evts")*makeList("HMS_PS"))/(makeList("charge")*makeList("curr_corr")*yield_dict["CPULT_phys"]*makeList("HMS_track"))

        
        
    yield_dict.update({"yield_HMS_scaler" : yield_HMS_scaler})
    yield_dict.update({"yield_HMS_notrack" : yield_HMS_notrack})
    yield_dict.update({"yield_HMS_track" : yield_HMS_track})
    yield_dict.update({"yield_HMS_CPULT_notrack" : yield_HMS_CPULT_notrack})
    yield_dict.update({"yield_HMS_CPULT_track" : yield_HMS_CPULT_track})

    for i,curr in enumerate(yield_dict["current"]):
        if abs((TLT_EDTM)[i] - (TLT_ELT)[i]) > ((abs(uncern_TLT_ELT)[i]) + abs((uncern_TLT)[i])) and TLT_ELT[i] > TLT_EDTM[i] :
            
            yield_dict["yield_HMS_notrack"][i] = (makeList("h_int_etotnorm_evts")[i]*makeList("HMS_PS")[i])/(makeList("charge")[i]*makeList("curr_corr")[i]*yield_dict["TLT_ELT"][i])
            yield_dict["yield_HMS_track"][i] = (makeList("h_int_etottracknorm_evts")[i]*makeList("HMS_PS")[i])/(makeList("charge")[i]*makeList("curr_corr")[i]*yield_dict["TLT_ELT"][i]*makeList("HMS_track")[i])
        else: pass 
    
    
    # define relative yield to the error weighted average of the absolute yield as opposed to the minimum current
    
    tot_hms_scaler = 0  #the total of each yield type, weighted by the reciprocal of their respective error
    tot_hms_tr = 0
    tot_hms_ntr = 0
    tot_hms_cpu_tr = 0
    tot_hms_cpu_ntr = 0
    
    tweight_hms_scaler = 0
    tweight_hms_tr = 0
    tweight_hms_ntr = 0
    tweight_hms_cpu_tr = 0
    tweight_hms_cpu_ntr = 0
    
    
    
    if makeList("HMS_PS").any() >0.5:
        uncern_yield_HMS_scaler = np.sqrt((yield_dict["yield_HMS_scaler"]**2)*((yield_dict["uncern_HMS_evts_scaler"]**2)+(yield_dict["uncern_CPULT_phys"]**2)+(yield_dict["uncern_charge"]/makeList("charge"))**2+(yield_dict["uncern_ELT2_HMS"])**2 + yield_dict["uncer_boilingCorr"] + yield_dict["uncer_rateHMSCorr"]))
        uncern_yield_HMS_track = np.sqrt(abs(yield_dict["yield_HMS_track"]**2)*((yield_dict["uncern_HMS_evts_scaler"]**2)+(yield_dict["uncern_CPULT_phys"]**2)+(yield_dict["uncern_charge"]/makeList("charge"))**2+ yield_dict["uncer_boilingCorr"] + yield_dict["uncer_rateHMSCorr"]))
        uncern_yield_HMS_notrack = np.sqrt(abs(yield_dict["yield_HMS_notrack"]**2)*((yield_dict["uncern_HMS_evts_scaler"]**2)+(yield_dict["uncern_CPULT_phys"]**2)+(yield_dict["uncern_charge"]/makeList("charge"))**2+ yield_dict["uncer_boilingCorr"] + yield_dict["uncer_rateHMSCorr"]))
        uncern_yield_HMS_CPULT_track = np.sqrt((yield_dict["yield_HMS_CPULT_track"]**2)*((yield_dict["uncern_HMS_evts_scaler"]**2)+(yield_dict["uncern_CPULT_phys"]**2)+(yield_dict["uncern_charge"]/makeList("charge"))**2))
        uncern_yield_HMS_CPULT_notrack = np.sqrt((yield_dict["yield_HMS_CPULT_notrack"]**2)*((yield_dict["uncern_HMS_evts_scaler"]**2)+(yield_dict["uncern_CPULT_phys"]**2)+(yield_dict["uncern_charge"]/makeList("charge"))**2))
    else:
        uncern_yield_HMS_scaler = .000000001
        uncern_yield_HMS_track = .00000001
        uncern_yield_HMS_notrack = .00000001
        uncern_yield_HMS_CPULT_track = .000000001
        uncern_yield_HMS_CPULT_notrack = .000000001
        
    if makeList("SHMS_PS").any() >0.5:
        uncern_yield_SHMS_scaler = np.sqrt((yield_dict["yield_SHMS_scaler"]**2)*((yield_dict["uncern_SHMS_evts_scaler"]**2)+(yield_dict["uncern_CPULT_phys"]**2)+(yield_dict["uncern_charge"]/makeList("charge"))**2 + yield_dict["uncer_boilingCorr"] + yield_dict["uncer_rateSHMSCorr"]))
        uncern_yield_SHMS_track = np.sqrt(abs(yield_dict["yield_SHMS_track"]**2)*((yield_dict["uncern_SHMS_evts_scaler"]**2)+(yield_dict["uncern_CPULT_phys"]**2)+(yield_dict["uncern_charge"]/makeList("charge"))**2+ yield_dict["uncer_boilingCorr"] + yield_dict["uncer_rateSHMSCorr"]))
        uncern_yield_SHMS_notrack = np.sqrt(abs(yield_dict["yield_SHMS_notrack"]**2)*((yield_dict["uncern_SHMS_evts_scaler"]**2)+(yield_dict["uncern_CPULT_phys"]**2)+(yield_dict["uncern_charge"]/makeList("charge"))**2+ yield_dict["uncer_boilingCorr"] + yield_dict["uncer_rateSHMSCorr"]))
        uncern_yield_SHMS_CPULT_track = np.sqrt((yield_dict["yield_SHMS_CPULT_track"]**2)*((yield_dict["uncern_SHMS_evts_scaler"]**2)+(yield_dict["uncern_CPULT_phys"]**2)+(yield_dict["uncern_charge"]/makeList("charge"))**2))
        uncern_yield_SHMS_CPULT_notrack = np.sqrt((yield_dict["yield_SHMS_CPULT_notrack"]**2)*((yield_dict["uncern_SHMS_evts_scaler"]**2)+(yield_dict["uncern_CPULT_phys"]**2)+(yield_dict["uncern_charge"]/makeList("charge"))**2))
    else:
        uncern_yield_SHMS_scaler = .000000001
        uncern_yield_SHMS_track = .00000001
        uncern_yield_SHMS_notrack = .00000001
        uncern_yield_SHMS_CPULT_track = .000000001
        uncern_yield_SHMS_CPULT_notrack = .000000001
            
    # Define relative yield relative to minimum current
    curr_tmp_shms = 0
    curr_tmp_hms = 0
    min_curr_shms = -1
    min_curr_hms = -1
    for i,curr in enumerate(yield_dict["current"]):
        if makeList("tot_events")[i] >= 10000:
            if makeList("SHMS_PS")[i] > 0.0:
                if len(yield_dict["current"]) <= 1:
                    min_curr_shms = yield_dict["current"][i]
                    break
                if curr_tmp_shms >= curr or curr_tmp_shms == 0:
                    min_curr_shms = yield_dict["current"][i]
                    curr_tmp_shms = curr
            else: # edge case for if SHMS not used in run
                min_curr_shms = 1
                min_yield_SHMS_scaler = 1
                min_yield_SHMS_notrack = 1
                min_yield_SHMS_track = 1
                min_yield_SHMS_CPULT_notrack = 1
                min_yield_SHMS_CPULT_track = 1
            if makeList("HMS_PS")[i] > 0.0:
                if len(yield_dict["current"]) <= 1:
                    min_curr_hms = yield_dict["current"][i]
                    break
                if curr_tmp_hms >= curr or curr_tmp_hms == 0:
                    min_curr_hms = yield_dict["current"][i]
                    curr_tmp_hms = curr
            else: # edge case for if HMS not used in run
                min_curr_hms = 1
                min_yield_HMS_scaler = 1
                min_yield_HMS_notrack = 1
                min_yield_HMS_track = 1
                min_yield_HMS_CPULT_notrack = 1
                min_yield_HMS_CPULT_track = 1
                
    for i,curr in enumerate(yield_dict["current"]):
        if (curr == min_curr_shms or curr == min_curr_hms):
            if makeList("SHMS_PS")[i] > 0.0:
                min_yield_SHMS_scaler = yield_dict["yield_SHMS_scaler"][i]
                min_yield_SHMS_notrack = yield_dict["yield_SHMS_notrack"][i]
                min_yield_SHMS_track = yield_dict["yield_SHMS_track"][i]
            if makeList("HMS_PS")[i] > 0.0:
                min_yield_HMS_scaler = yield_dict["yield_HMS_scaler"][i]
                min_yield_HMS_notrack = yield_dict["yield_HMS_notrack"][i]
                min_yield_HMS_track = yield_dict["yield_HMS_track"][i]
            #print("run a min current")
            min_yield_COIN_scaler = yield_dict["yield_COIN_scaler"][i]
            min_yield_COIN_notrack = yield_dict["yield_COIN_notrack"][i]
            min_yield_COIN_track = yield_dict["yield_COIN_track"][i]
    yield_dict.update({"min_yield_HMS_scaler" : min_yield_HMS_scaler})                
    yield_dict.update({"min_yield_HMS_notrack" : min_yield_HMS_notrack})
    yield_dict.update({"min_yield_HMS_track" : min_yield_HMS_track})
    yield_dict.update({"min_yield_SHMS_scaler" : min_yield_SHMS_scaler})                
    yield_dict.update({"min_yield_SHMS_notrack" : min_yield_SHMS_notrack})
    yield_dict.update({"min_yield_SHMS_track" : min_yield_SHMS_track})
    yield_dict.update({"min_yield_COIN_scaler": min_yield_COIN_scaler})
    yield_dict.update({"min_yield_COIN_notrack": min_yield_COIN_notrack})
    yield_dict.update({"min_yield_COIN_track": min_yield_COIN_track})

    # Define relative yield relative to minimum current
    for i,curr in enumerate(yield_dict["current"]):
        if curr ==  min_curr_shms:
            min_yield_SHMS_CPULT_notrack = yield_dict["yield_SHMS_CPULT_notrack"][i]
            min_yield_SHMS_CPULT_track = yield_dict["yield_SHMS_CPULT_track"][i]             
        if curr ==  min_curr_hms:
            min_yield_HMS_CPULT_notrack = yield_dict["yield_HMS_CPULT_notrack"][i]
            min_yield_HMS_CPULT_track = yield_dict["yield_HMS_CPULT_track"][i]             
    yield_dict.update({"min_yield_HMS_CPULT_notrack" : min_yield_HMS_CPULT_notrack})
    yield_dict.update({"min_yield_HMS_CPULT_track" : min_yield_HMS_CPULT_track})
    yield_dict.update({"min_yield_SHMS_CPULT_notrack" : min_yield_SHMS_CPULT_notrack})
    yield_dict.update({"min_yield_SHMS_CPULT_track" : min_yield_SHMS_CPULT_track})

    if makeList("SHMS_PS").any() != 0:
        yieldRel_SHMS_scaler = yield_dict["yield_SHMS_scaler"]/yield_dict["min_yield_SHMS_scaler"]
        yieldRel_SHMS_notrack = yield_dict["yield_SHMS_notrack"]/yield_dict["min_yield_SHMS_notrack"]
        yieldRel_SHMS_track = yield_dict["yield_SHMS_track"]/yield_dict["min_yield_SHMS_track"]
        yield_dict.update({"yieldRel_SHMS_scaler" : yieldRel_SHMS_scaler})
        yield_dict.update({"yieldRel_SHMS_notrack" : yieldRel_SHMS_notrack})
        yield_dict.update({"yieldRel_SHMS_track" : yieldRel_SHMS_track})

        yieldRel_SHMS_CPULT_notrack = yield_dict["yield_SHMS_CPULT_notrack"]/yield_dict["min_yield_SHMS_CPULT_notrack"]
        yieldRel_SHMS_CPULT_track = yield_dict["yield_SHMS_CPULT_track"]/yield_dict["min_yield_SHMS_CPULT_track"]
        yield_dict.update({"yieldRel_SHMS_CPULT_notrack" : yieldRel_SHMS_CPULT_notrack})
        yield_dict.update({"yieldRel_SHMS_CPULT_track" : yieldRel_SHMS_CPULT_track})
            
        uncern_yieldRel_SHMS_scaler = np.sqrt((yield_dict["uncern_SHMS_evts_scaler"]**2)+(yield_dict["uncern_CPULT_phys"]**2)+(yield_dict["uncern_charge"]/makeList("charge"))**2 + yield_dict["uncer_boilingCorr"] + yield_dict["uncer_rateSHMSCorr"])
        uncern_yieldRel_SHMS_notrack = np.sqrt((yield_dict["uncern_SHMS_evts_notrack"]**2)+(yield_dict["uncern_TLT_ELT"])+(yield_dict["uncern_charge"]/makeList("charge"))**2 + yield_dict["uncer_boilingCorr"] + yield_dict["uncer_rateSHMSCorr"])
        uncern_yieldRel_SHMS_track =  np.sqrt((yield_dict["uncern_SHMS_evts_notrack"]**2)+(yield_dict["uncern_TLT_ELT"])+(makeList("SHMS_track_uncern")**2)+((yield_dict["uncern_charge"]/makeList("charge"))**2) + yield_dict["uncer_boilingCorr"] + yield_dict["uncer_rateSHMSCorr"])
        yield_dict.update({"uncern_yieldRel_SHMS_scaler" : uncern_yieldRel_SHMS_scaler})
        yield_dict.update({"uncern_yieldRel_SHMS_notrack" : uncern_yieldRel_SHMS_notrack})
        yield_dict.update({"uncern_yieldRel_SHMS_track" : uncern_yieldRel_SHMS_track})
            
        uncern_yieldRel_SHMS_CPULT_notrack = (yield_dict["uncern_SHMS_evts_notrack"]+yield_dict["uncern_CPULT_phys"]+yield_dict["uncern_charge"]/makeList("charge"))
        uncern_yieldRel_SHMS_CPULT_track =  (yield_dict["uncern_SHMS_evts_notrack"]+yield_dict["uncern_CPULT_phys"]+makeList("SHMS_track_uncern")+yield_dict["uncern_charge"]/makeList("charge"))
        yield_dict.update({"uncern_yieldRel_SHMS_CPULT_notrack" : uncern_yieldRel_SHMS_CPULT_notrack})
        yield_dict.update({"uncern_yieldRel_SHMS_CPULT_track" : uncern_yieldRel_SHMS_CPULT_track})
            
    else: 
        yield_dict.update({"yieldRel_SHMS_scaler" :0})
        yield_dict.update({"yieldRel_SHMS_notrack" : 0})
        yield_dict.update({"yieldRel_SHMS_track" : 0})
        yield_dict.update({"yieldRel_SHMS_CPULT_notrack" : 0})
        yield_dict.update({"yieldRel_SHMS_CPULT_track" : 0})
            
        yield_dict.update({"uncern_yieldRel_SHMS_scaler" : 0})
        yield_dict.update({"uncern_yieldRel_SHMS_notrack" : 0})
        yield_dict.update({"uncern_yieldRel_SHMS_track" : 0})
        
        yield_dict.update({"uncern_yieldRel_SHMS_CPULT_notrack" : 0})
        yield_dict.update({"uncern_yieldRel_SHMS_CPULT_track" : 0})
            
    if makeList("HMS_PS").any() != 0:
        yieldRel_HMS_scaler = yield_dict["yield_HMS_scaler"]/yield_dict["min_yield_HMS_scaler"]
        yieldRel_HMS_notrack = yield_dict["yield_HMS_notrack"]/yield_dict["min_yield_HMS_notrack"]
        yieldRel_HMS_track = yield_dict["yield_HMS_track"]/yield_dict["min_yield_HMS_track"]
        yield_dict.update({"yieldRel_HMS_scaler" : yieldRel_HMS_scaler})
        yield_dict.update({"yieldRel_HMS_notrack" : yieldRel_HMS_notrack})
        yield_dict.update({"yieldRel_HMS_track" : yieldRel_HMS_track})

        yieldRel_HMS_CPULT_notrack = yield_dict["yield_HMS_CPULT_notrack"]/yield_dict["min_yield_HMS_CPULT_notrack"]
        yieldRel_HMS_CPULT_track = yield_dict["yield_HMS_CPULT_track"]/yield_dict["min_yield_HMS_CPULT_track"]
        yield_dict.update({"yieldRel_HMS_CPULT_notrack" : yieldRel_HMS_CPULT_notrack})
        yield_dict.update({"yieldRel_HMS_CPULT_track" : yieldRel_HMS_CPULT_track})
            
        uncern_yieldRel_HMS_scaler = np.sqrt((yield_dict["uncern_HMS_evts_scaler"]**2)+(yield_dict["uncern_charge"]/makeList("charge"))**2 +(yield_dict["uncern_ELT2_HMS"])**2 + yield_dict["uncer_boilingCorr"]  + yield_dict["uncer_rateHMSCorr"])
        uncern_yieldRel_HMS_notrack = np.sqrt((yield_dict["uncern_HMS_evts_notrack"]**2)+(yield_dict["uncern_TLT_ELT"])+(yield_dict["uncern_charge"]/makeList("charge"))**2 + yield_dict["uncer_boilingCorr"]  + yield_dict["uncer_rateHMSCorr"])
        uncern_yieldRel_HMS_track =  np.sqrt((yield_dict["uncern_HMS_evts_notrack"]**2)+(yield_dict["uncern_TLT_ELT"])+(makeList("HMS_track_uncern")**2)+(yield_dict["uncern_charge"]/makeList("charge"))**2 + yield_dict["uncer_boilingCorr"]  + yield_dict["uncer_rateHMSCorr"])
        yield_dict.update({"uncern_yieldRel_HMS_scaler" : uncern_yieldRel_HMS_scaler})
        yield_dict.update({"uncern_yieldRel_HMS_notrack" : uncern_yieldRel_HMS_notrack})
        yield_dict.update({"uncern_yieldRel_HMS_track" : uncern_yieldRel_HMS_track})
            
        uncern_yieldRel_HMS_CPULT_notrack = (yield_dict["uncern_HMS_evts_notrack"]+yield_dict["uncern_CPULT_phys"])
        uncern_yieldRel_HMS_CPULT_track =  (yield_dict["uncern_HMS_evts_notrack"]+yield_dict["uncern_CPULT_phys"]+makeList("HMS_track_uncern"))
        yield_dict.update({"uncern_yieldRel_HMS_CPULT_notrack" : uncern_yieldRel_HMS_CPULT_notrack})
        yield_dict.update({"uncern_yieldRel_HMS_CPULT_track" : uncern_yieldRel_HMS_CPULT_track})
    else:
        yield_dict.update({"yieldRel_HMS_scaler" : 0})
        yield_dict.update({"yieldRel_HMS_notrack" : 0})
        yield_dict.update({"yieldRel_HMS_track" : 0})
        yield_dict.update({"yieldRel_HMS_CPULT_notrack" : 0})
        yield_dict.update({"yieldRel_HMS_CPULT_track" : 0})
        
        yield_dict.update({"uncern_yieldRel_HMS_scaler" : 0})
        yield_dict.update({"uncern_yieldRel_HMS_notrack" : 0})
        yield_dict.update({"uncern_yieldRel_HMS_track" : 0})
            
        yield_dict.update({"uncern_yieldRel_HMS_CPULT_notrack" : 0})
        yield_dict.update({"uncern_yieldRel_HMS_CPULT_track" : 0})
        
    
    yieldRel_COIN_scaler = yield_dict["yield_COIN_scaler"]/yield_dict["min_yield_COIN_scaler"]
    yieldRel_COIN_notrack = yield_dict["yield_COIN_notrack"]/yield_dict["min_yield_COIN_notrack"]
    yieldRel_COIN_track = yield_dict["yield_COIN_track"]/yield_dict["min_yield_COIN_track"]
    
    print("\n\n******************")
    print("yieldRel_COIN_track:")
    print(yieldRel_COIN_track)
    print("\n\n*******************")
    yield_dict.update({"yieldRel_COIN_scaler": yieldRel_COIN_scaler})
    yield_dict.update({"yieldRel_COIN_notrack": yieldRel_COIN_notrack})
    yield_dict.update({"yieldRel_COIN_track": yieldRel_COIN_track})
        
    #for i,psmod in enumerate(yield_dict["PS_mod"]):
    #    if psmod > 0:
    #        makeList("sent_edtm_PS")[i] = yield_dict["sent_edtm_PS_corr"][i]
        #else:
         #   sent_edtm_PS_final[i] = makeList("sent_edtm_PS")[i]            
    #yield_dict.update({"sent_edtm_PS_final" : sent_edtm_PS_final})

    
    # Restructure dictionary to dataframe format so it matches lumi_data
    yield_table = pd.DataFrame(yield_dict, columns=yield_dict.keys())
    yield_table = yield_table.reindex(sorted(yield_table.columns), axis=1)

    return yield_table

def mergeDicts():
    '''
    Merge dictionaries/dataframes, convert to dataframe and sort
    '''
    yield_data = calc_yield()
    # data = {**lumi_data, **yield_data} # only python 3.5+
    
    for key, val in lumi_data.items():
        lumi_data[key] = val

    datadict = {}
    for d in (lumi_data, yield_data): 
        datadict.update(d)
    data = {i : datadict[i] for i in sorted(datadict.keys())}

    #data = {x:y for x,y in data.items() if print("HERE",y!=0)}
    #data = {x:y.replace(0,np.nan,inplace=True) for x,y in data.items()}

    table  = pd.DataFrame(data, columns=data.keys())
    table = table.reindex(sorted(table.columns), axis=1)

    # Replace zeros with NaN
    table.replace(0,np.nan,inplace=True)

    return table

################################################################################################################################################

def plot_yield():
    '''
    Plot yields and various other analysis plots
    '''
    
    # determine which scaler is being used by pulling out the prescales
    prescales = {
        "PS1" : makeList("PS1"),
        "PS2" : makeList("PS2"),
        "PS3" : makeList("PS3"),
        "PS4" : makeList("PS4"),
        "PS5" : makeList("PS5"),
        "PS6" : makeList("PS6"),
    }
    
    HMSscaler = "OFF"
    SHMSscaler = "OFF"
    COINscaler = "OFF"
    
    if (not (prescales['PS1'][1] == None)):
         SHMSscaler = "3/4 (PS1)"
    if (not (prescales['PS2'][1] == None)):
         SHMSscaler = "ElReal (PS2)"
    if (not (prescales['PS3'][1] == None)):
         HMSscaler = "3/4 (PS3)"
    if (not (prescales['PS4'][1] == None)):
         HMSscaler = "ElReal (PS4)"
    if (not (prescales['PS5'][1] == None)):
         COINscaler = "3/4xElReal (PS5)"
    if (not (prescales['PS6'][1] == None)):
         COINscaler = "3/4x3/4 (PS6)"
    
    yield_data = mergeDicts()

    # Remove runs with bad TLT, print statement only
    #for i,row in yield_data.iterrows():
    #    if row['TLT'] <= 0.75 or row['TLT'] > 1.02:
    #        '''
    #        print(''
    #        Removing {0:.0f} because TLT is low...
    #        TLT={1:0.2f} %
    #        CPULT={2:0.2f} %
    #        ''.format(row["run number"],row["TLT"]*100,row["CPULT_phys"]*100))
    #        '''
    #        print('''
    #        Bad TLT...
    #        TLT={0:0.2f} %
    #        CPULT={1:0.2f} %
    #        '''.format(row["TLT"]*100,row["CPULT_phys"]*100))
    #
    #    if row['time'] < 60.0:
    #        print('''
    #        Removing {0:.0f} because beam on time is short...
    #        time={1:0.1f} s
    #        '''.format(row["run number"],row["time"]))

    # Remove runs with bad TLT or short on beam time
    #yield_data = yield_data[ (yield_data['TLT'] >= 0.75) & (yield_data['CPULT_phys'] >= 0.75) ].reset_index(drop=True)
    #yield_data = yield_data[ (yield_data['TLT'] < 1.02) & (yield_data['CPULT_phys] < 1.02)].reset_index(drop=True)
    #yield_data = yield_data[ (yield_data['time'] >= 60.0) ].reset_index(drop=True)

    for i, val in enumerate(yield_data["run number"]):
        print("Run numbers:",yield_data["run number"][i],"Current Values:",yield_data["current"][i])

    
    #def fit_func(x, m, b):
       # y = (m/b)*x + 1.0
       # return y
    '''
    def linear_plot(x, y, xerr=None, yerr=None, xvalmax=100):
        # Remove NaN values from list
        c_arr = [[nx,ny] for nx,ny in zip(x,y) if str(ny) != 'nan']
        x_c = np.array([i[0] for i in c_arr])
        y_c = np.array([i[1] for i in c_arr])

        # Define slope and intercept
        if len(x_c) < 4:
            popt, _ = curve_fit(fit_func, x_c, y_c)
            m, b = popt
            m_err = np.std(y_c - fit_func(x_c, m, b)) / np.sqrt(len(x_c))
            b_err = np.std(y_c - fit_func(x_c, m, b)) / np.sqrt(len(x_c))
        else:
            popt, pcov = np.polyfit(x_c, y_c, deg=1, cov=True)
            m, b = popt
            m_err, b_err = np.sqrt(np.diag(pcov))

        # Find chi-squared
        res = y_c - fit_func(x_c, m, b)
        chisq = np.sum((res/(y_c*yerr))**2) if yerr is not None else np.sum(res**2)

        # Plot fit from axis
        x_fit = np.linspace(0,xvalmax)
        y_fit = fit_func(x_fit, m/b, 1.0)
        plt.plot(x_fit, y_fit, color='green', label='{0}={1:0.2e}*{2}+1.00\n{3}={4:0.2e}\n{5}={6:0.2e}\n{7}={8:0.2e}'.format(r'Y/$Y_0$',m/b,r'$I_b$',r'$\chi^2$',chisq,r'$m_0$',(m/b),r'$\delta m_0$',(m_err/b_err)), zorder=5)

        return m/b
    '''

    '''def linear_plot(x, y, xerr=None, yerr=None, xvalmax=100):
        # Remove NaN values from list
        c_arr = [[nx,ny] for nx,ny in zip(x,y) if str(ny) != 'nan']
        x_c = np.array([i[0] for i in c_arr])
        y_c = np.array([i[1] for i in c_arr])

        # perform weighted least squares regression
        results = sm.WLS(y_c[:, np.newaxis], sm.add_constant(x_c[:, np.newaxis]), weights=1.0/yerr[:, np.newaxis]**2).fit()

        m = results.params[1]
        b = results.params[0]
        
        # Find chi-squared
        res = y_c - results.predict(sm.add_constant(x_c))
        chisq = np.sum((res/(y_c*yerr))**2) if yerr is not None else np.sum(res**2)

        # Plot fit from axis
        x_fit = np.linspace(0,xvalmax)
        y_fit = results.predict(sm.add_constant(x_fit))
        plt.plot(x_fit, y_fit, color='green', label='{0}={1:0.2e}*{2}+1.00\n{3}={4:0.2e}\n{5}={6:0.2e}'.format(r'Y/$Y_0$',m/b,r'$I_b$',r'$\chi^2$',chisq,r'$m_0$',(m/b)), zorder=5)

        return m/b  '''  

    def linfunc(x, int, slope):
        y = int + slope*x
        return y
    
    def linfunc2(x,slope):
    
        y = slope*x + 1.0
        return y
    
    def constfunc(x,slope):
        
        y = slope
        return y
    
    def quadfunc(x,int,linear,quadratic):
        
        y = int  + (linear)*(x)+ (quadratic/6400)*((x)**(2))
        return y
    
   # def expfunc(x, b, a):
    #    y = b^(a*x)
      #  return y
    
        
    
    
    #########################################################################################################################################################

    relYieldPlot = plt.figure(figsize=(12,8))
    if(not NO_HMS):
        #HMS plot scaler
        plt.subplot(3,3,1)    
        plt.grid(zorder=1)
        plt.xlim(0,100)
        plt.ylim(0.9,1.1)
        plt.plot([0,100], [1,1], 'r-',zorder=2)
        
        
        if(np.isnan(yield_data["yieldRel_HMS_scaler"][0])== False):
            a_fit_HMS_scalerVScurrent, cov_HMS_scalerVScurrent = curve_fit(linfunc,yield_data["current"], yield_data["yieldRel_HMS_scaler"], sigma=yield_data["uncern_yieldRel_HMS_scaler"], absolute_sigma = True)
            inter_HMS_scalerVScurrent = a_fit_HMS_scalerVScurrent[0]
            slope_HMS_scalerVScurrent = a_fit_HMS_scalerVScurrent[1]
            d_inter_HMS_scalerVScurrent = np.sqrt(cov_HMS_scalerVScurrent[0][0])
            d_slope_HMS_scalerVScurrent = np.sqrt(cov_HMS_scalerVScurrent[1][1])
            
            a_fit_HMS_scalerVScurrent2, cov_HMS_scalerVScurrent2 = curve_fit(linfunc2,yield_data["current"], yield_data["yieldRel_HMS_scaler"], sigma=yield_data["uncern_yieldRel_HMS_scaler"], absolute_sigma = True)
           
            slope_HMS_scalerVScurrent2 = a_fit_HMS_scalerVScurrent2[0]
        
            d_slope_HMS_scalerVScurrent2= np.sqrt(cov_HMS_scalerVScurrent2[0][0])
        
        plt.errorbar(yield_data["current"],yield_data["yieldRel_HMS_scaler"],yerr=yield_data["yieldRel_HMS_scaler"]*yield_data["uncern_yieldRel_HMS_scaler"],color='black',linestyle='None',zorder=3,label="_nolegend_")
        plt.scatter(yield_data["current"],yield_data["yieldRel_HMS_scaler"],color='blue',zorder=4,label="_nolegend_")
        #    yield_data["m0_curr_HMS_scaler"] = linear_plot(yield_data["current"],yield_data["yieldRel_HMS_scaler"],None,yield_data["uncern_yieldRel_HMS_scaler"])
        
       
        
        if(np.isnan(yield_data["yieldRel_HMS_scaler"][0])== False):
            yfit = (slope_HMS_scalerVScurrent)*yield_data["current"] + inter_HMS_scalerVScurrent
            yfit2 = (slope_HMS_scalerVScurrent2)*yield_data["current"] + 1.0
            plt.plot(yield_data["current"], yfit, color = 'green', label = "slope = %f +/- %f" %(slope_HMS_scalerVScurrent, d_slope_HMS_scalerVScurrent) + "\n intercept = %f +/- %f" %(inter_HMS_scalerVScurrent, d_inter_HMS_scalerVScurrent))

        
        plt.ylabel('Rel. Yield %s' % (str(HMSscaler)), fontsize=16)
        plt.xlabel('Current [uA]', fontsize =16)
        plt.legend(prop={'size' :7})
        if target == 'LD2' :
            plt.title('HMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
        elif target == 'LH2' :
            plt.title('HMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
        else :
            plt.title('HMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    
        #HMS plot no track
        plt.subplot(3,3,2)    
        plt.grid(zorder=1)
        plt.xlim(0,100)
        plt.ylim(0.9,1.1)
        plt.plot([0,100], [1,1], 'r-',zorder=2)
        
        plt.errorbar(yield_data["current"], yield_data["yieldRel_HMS_notrack"], yerr=yield_data["yieldRel_HMS_notrack"]*yield_data["uncern_yieldRel_HMS_notrack"], color='black', linestyle='None', zorder=3, label="_nolegend_")
        plt.scatter(yield_data["current"],yield_data["yieldRel_HMS_notrack"],color='blue',zorder=4,label="_nolegend_")
        #    yield_data["m0_curr_HMS_notrack"] = linear_plot(yield_data["current"],yield_data["yieldRel_HMS_CPULT_notrack"],None,yield_data["uncern_yieldRel_HMS_CPULT_notrack"])      
        plt.ylabel('Rel. Yield no track', fontsize=16)
        plt.xlabel('Current [uA]', fontsize =16)
        plt.legend(prop={'size' :7})
        if target == 'LD2' :
            plt.title('HMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
        elif target == 'LH2' :
            plt.title('HMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
        else :
            plt.title('HMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
        
        #HMS plot track
        plt.subplot(3,3,3)    
        plt.grid(zorder=1)
        plt.xlim(0,100)
        plt.ylim(0.9,1.1)
        plt.plot([0,100], [1,1], 'r-',zorder=2)
        plt.errorbar(yield_data["current"],yield_data["yieldRel_HMS_track"],yerr=yield_data["yieldRel_HMS_track"]*yield_data["uncern_yieldRel_HMS_track"],color='black',linestyle='None',zorder=3,label="_nolegend_")
        plt.scatter(yield_data["current"],yield_data["yieldRel_HMS_track"],color='blue',zorder=4,label="_nolegend_")
        #    yield_data["m0_curr_HMS_track"] = linear_plot(yield_data["current"],yield_data["yieldRel_HMS_track"],None,yield_data["uncern_yieldRel_HMS_track"])
        #plt.errorbar(yield_data["current"],yield_data["yieldRel_HMS_CPULT_track"],yerr=yield_data["yieldRel_HMS_CPULT_track"]*yield_data["uncern_yieldRel_HMS_CPULT_track"],color='black',linestyle='None',zorder=5)
        #plt.scatter(yield_data["current"],yield_data["yieldRel_HMS_CPULT_track"],color='red',zorder=6)
        plt.ylabel('Rel. Yield track', fontsize=16)
        plt.xlabel('Current [uA]', fontsize =16)
        plt.legend(prop={'size' :7})
        if target == 'LD2' :
            plt.title('HMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
        elif target == 'LH2' :
            plt.title('HMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
        else :
            plt.title('HMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    
    if(not NO_SHMS):    
        #SHMS plot scaler
        plt.subplot(3,3,4)    
        plt.grid(zorder=1)
        plt.xlim(0,100)
        plt.ylim(0.9,1.1)
        plt.plot([0,100], [1,1], 'r-',zorder=2)
        
        plt.errorbar(yield_data["current"],yield_data["yieldRel_SHMS_scaler"],yerr=yield_data["yieldRel_SHMS_scaler"]*yield_data["uncern_yieldRel_SHMS_scaler"],color='black',linestyle='None',zorder=3,label="_nolegend_")
        plt.scatter(yield_data["current"],yield_data["yieldRel_SHMS_scaler"],color='blue',zorder=4,label="_nolegend_")
        #    yield_data["m0_curr_SHMS_scaler"] = linear_plot(yield_data["current"],yield_data["yieldRel_SHMS_scaler"],None,yield_data["uncern_yieldRel_SHMS_scaler"])

        plt.ylabel('Rel. Yield %s' % (str(SHMSscaler)), fontsize=16)
        plt.xlabel('Current [uA]', fontsize =16)
        plt.legend(prop={'size' :7})
        if target == 'LD2' :
            plt.title('SHMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
        elif target == 'LH2' :
            plt.title('SHMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
        else :
            plt.title('SHMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    
        #SHMS plot no track
        plt.subplot(3,3,5)    
        plt.grid(zorder=1)
        plt.xlim(0,100)
        plt.ylim(0.9,1.1)
        plt.plot([0,100], [1,1], 'r-',zorder=2)
        
        plt.errorbar(yield_data["current"],yield_data["yieldRel_SHMS_notrack"],yerr=yield_data["yieldRel_SHMS_notrack"]*yield_data["uncern_yieldRel_SHMS_notrack"],color='black',linestyle='None',zorder=3,label="_nolegend_")
        plt.scatter(yield_data["current"],yield_data["yieldRel_SHMS_notrack"],color='blue',zorder=4,label="_nolegend_")
        #    yield_data["m0_curr_SHMS_notrack"] = linear_plot(yield_data["current"],yield_data["yieldRel_SHMS_notrack"],None,yield_data["uncern_yieldRel_SHMS_notrack"])
        #plt.errorbar(yield_data["current"],yield_data["yieldRel_SHMS_CPULT_notrack"],yerr=yield_data["yieldRel_SHMS_CPULT_notrack"]*yield_data["uncern_yieldRel_SHMS_CPULT_notrack"],color='black',linestyle='None',zorder=5)
        #plt.scatter(yield_data["current"],yield_data["yieldRel_SHMS_CPULT_notrack"],color='red',zorder=6)
        
        plt.ylabel('Rel. Yield no track', fontsize=16)
        plt.xlabel('Current [uA]', fontsize =16)
        plt.legend(prop={'size' :7})
        if target == 'LD2' :
            plt.title('SHMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
        elif target == 'LH2' :
            plt.title('SHMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
        else :
            plt.title('SHMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
            
        #SHMS plot track
        plt.subplot(3,3,6)    
        plt.grid(zorder=1)
        plt.xlim(0,100)
        plt.ylim(0.9,1.1)
        plt.plot([0,100], [1,1], 'r-',zorder=2)
        
        plt.errorbar(yield_data["current"],yield_data["yieldRel_SHMS_track"],yerr=yield_data["yieldRel_SHMS_track"]*yield_data["uncern_yieldRel_SHMS_track"],color='black',linestyle='None',zorder=3,label="_nolegend_")
        plt.scatter(yield_data["current"],yield_data["yieldRel_SHMS_track"],color='blue',zorder=4,label="_nolegend_")
        #    yield_data["m0_curr_SHMS_track"] = linear_plot(yield_data["current"],yield_data["yieldRel_SHMS_track"]0,None,yield_data["uncern_yieldRel_SHMS_track"])
        #plt.errorbar(yield_data["current"],yield_data["yieldRel_SHMS_CPULT_track"],yerr=yield_data["yieldRel_SHMS_CPULT_track"]*yield_data["uncern_yieldRel_SHMS_CPULT_track"],color='black',linestyle='None',zorder=5)
        #plt.scatter(yield_data["current"],yield_data["yieldRel_SHMS_CPULT_track"],color='red',zorder=6)
        
        plt.ylabel('Rel. Yield track', fontsize=16)
        plt.xlabel('Current [uA]', fontsize =16)
        plt.legend(prop={'size' :7})
        if target == 'LD2' :
            plt.title('SHMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
        elif target == 'LH2' :
            plt.title('SHMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
        else :
            plt.title('SHMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    
    # Start Coin plots
    ## Scaler ##
    plt.subplot(3,3,7)    
    plt.grid(zorder=1)
    plt.xlim(0,100)
    plt.ylim(0.9,4)
    plt.plot([0,100], [1,1], 'r-',zorder=2)
    plt.errorbar(yield_data["current"],yield_data["yieldRel_COIN_scaler"],yerr=yield_data["yieldRel_COIN_scaler"]*yield_data["uncern_yieldRel_COIN_scaler"],color='black',linestyle='None',zorder=3,label="_nolegend_")
    plt.scatter(yield_data["current"],yield_data["yieldRel_COIN_scaler"],color='blue',zorder=4,label="_nolegend_")
    
    plt.ylabel('Rel. Yield %s' % (str(COINscaler)), fontsize=16)
    plt.xlabel('Current [uA]', fontsize =16)
    plt.legend(prop={'size' :7})    
    if target == 'LD2' :
        plt.title('COIN LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('COIN LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('COIN Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
            
    ## No Track ##        
    plt.subplot(3,3,8)    
    plt.grid(zorder=1)
    plt.xlim(0,100)
    plt.ylim(0.85,1.05)
    plt.plot([0,100], [1,1], 'r-',zorder=2)
    plt.errorbar(yield_data["current"],yield_data["yieldRel_COIN_notrack"],yerr=yield_data["yieldRel_COIN_notrack"]*yield_data["uncern_yieldRel_COIN_notrack"],color='black',linestyle='None',zorder=3,label="_nolegend_")
    plt.scatter(yield_data["current"],yield_data["yieldRel_COIN_notrack"],color='blue',zorder=4,label="_nolegend_")
    
    plt.ylabel('Rel. Yield no track', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =16)
    plt.legend(prop={'size' :7})    
    if target == 'LD2' :
        plt.title('COIN LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('COIN LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('COIN Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    
    ## Track ##
    plt.subplot(3,3,9)    
    plt.grid(zorder=1)
    plt.xlim(0,100)
    plt.ylim(0.85,1.05)
    plt.plot([0,100], [1,1], 'r-',zorder=2)
    plt.errorbar(yield_data["current"],yield_data["yieldRel_COIN_track"],yerr=yield_data["yieldRel_COIN_track"]*yield_data["uncern_yieldRel_COIN_track"],color='black',linestyle='None',zorder=3,label="_nolegend_")
    plt.scatter(yield_data["current"],yield_data["yieldRel_COIN_track"],color='blue',zorder=4,label="_nolegend_")
    
    plt.ylabel('Rel. Yield track', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =16)
    plt.legend(prop={'size' :7})    
    if target == 'LD2' :
        plt.title('COIN LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('COIN LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('COIN Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
                
    plt.tight_layout()
    plt.savefig(SCRIPTPATH+'/luminosity/OUTPUTS/plots/Yield_%s_%s.png' % (out_f.split("yield_data_")[1].replace(".csv",""),"relYieldPlot"))
        
    #########################################################################################################################################################
    ##              Rate Depedance Plots        #############################################################################################################
    #########################################################################################################################################################
    
    RatePlotFontSize = 8
    raterelYieldPlot = plt.figure(figsize=(12,8))

    if(not NO_HMS):
        #HMS plot scaler
        plt.subplot(2,3,1)    
        plt.grid(zorder=1)
        plt.xlim(0,(max(yield_data["rate_HMS"])/1000)+5)
        plt.ylim(0.9,1.1)
        plt.plot([0,(max(yield_data["rate_HMS"])/1000)+5], [1,1], 'r-',zorder=2)
        
        if(np.isnan(yield_data["yieldRel_HMS_scaler"][0])== False):
            a_fit_HMS_scalerVSrate, cov_HMS_scalerVSrate = curve_fit(linfunc,yield_data["rate_HMS"]/1000, yield_data["yieldRel_HMS_scaler"], sigma=yield_data["uncern_yieldRel_HMS_scaler"], absolute_sigma = True)
            inter_HMS_scalerVSrate = a_fit_HMS_scalerVSrate[0]
            slope_HMS_scalerVSrate = a_fit_HMS_scalerVSrate[1]
            d_inter_HMS_scalerVSrate = np.sqrt(cov_HMS_scalerVSrate[0][0])
            d_slope_HMS_scalerVSrate = np.sqrt(cov_HMS_scalerVSrate[1][1])
        
        plt.errorbar(yield_data["rate_HMS"]/1000,yield_data["yieldRel_HMS_scaler"],yerr=yield_data["yieldRel_HMS_scaler"]*yield_data["uncern_yieldRel_HMS_scaler"],color='black',linestyle='None',zorder=3,label="_nolegend_")
        plt.scatter(yield_data["rate_HMS"]/1000,yield_data["yieldRel_HMS_scaler"],color='blue',zorder=4,label="_nolegend_")
        #yield_data["m0_rate_HMS_scaler"] = linear_plot(yield_data["rate_HMS"]/1000,yield_data["yieldRel_HMS_scaler"],None,yield_data["uncern_yieldRel_HMS_scaler"],xvalmax=max((yield_data["rate_HMS"])/1000)+5)
        
        if(np.isnan(yield_data["yieldRel_HMS_scaler"][0])== False):
            yfit = (slope_HMS_scalerVSrate)*(yield_data["rate_HMS"]/1000) + inter_HMS_scalerVSrate
            plt.plot(yield_data["rate_HMS"]/1000, yfit, color = 'green', label = "slope = %f +/- %f" %(slope_HMS_scalerVSrate, d_slope_HMS_scalerVSrate) + "\n intercept = %f +/- %f" %(inter_HMS_scalerVSrate, d_inter_HMS_scalerVSrate))
    
        
        print("debug HMS rate ploting:")
        print(yield_data["rate_HMS"]/1000)
        print(yield_data["yieldRel_HMS_scaler"])
        
        #plot line of best fit
        #coeff_HMS_scalerVSrate = np.polyfit(yield_data["rate_HMS"], yield_data["yieldRel_HMS_scaler"], 1)
        #poly_HMS_scalerVSrate = np.poly1d(coeff_HMS_scalerVSrate)
        #plt.plot(yield_data["rate_HMS"]/1000, poly_HMS_scalerVSrate(yield_data["yieldRel_HMS_scaler"]), color='green', label='Line of Best Fit')
       
        plt.ylabel('Rel. Yield %s' % (str(HMSscaler)), fontsize=RatePlotFontSize)
        plt.xlabel('HMS %s Rate [kHz]' % (str(HMSscaler)), fontsize=RatePlotFontSize)
        plt.legend(prop={'size' :7})
        if target == 'LD2' :
            plt.title('HMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize=RatePlotFontSize)
        elif target == 'LH2' :
            plt.title('HMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize=RatePlotFontSize)
        else :
            plt.title('HMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize=RatePlotFontSize)
    
        #HMS plot no track
        plt.subplot(2,3,2)    
        plt.grid(zorder=1)
        plt.xlim(0,(max(yield_data["rate_HMS"])/1000)+5)
        plt.ylim(0.9,1.1)
        plt.plot([0,(max(yield_data["rate_HMS"])/1000)+5], [1,1], 'r-',zorder=2)
        
        if(np.isnan(yield_data["yieldRel_HMS_notrack"][0])== False):
            a_fit_HMS_ntrVSrate, cov_HMS_ntrVSrate = curve_fit(linfunc,yield_data["rate_HMS"]/1000, yield_data["yieldRel_HMS_notrack"], sigma=yield_data["uncern_yieldRel_HMS_notrack"], absolute_sigma = True)
            inter_HMS_ntrVSrate = a_fit_HMS_ntrVSrate[0]
            slope_HMS_ntrVSrate = a_fit_HMS_ntrVSrate[1]
            d_inter_HMS_ntrVSrate = np.sqrt(cov_HMS_ntrVSrate[0][0])
            d_slope_HMS_ntrVSrate = np.sqrt(cov_HMS_ntrVSrate[1][1])
        
        
        plt.errorbar(yield_data["rate_HMS"]/1000,yield_data["yieldRel_HMS_notrack"],yerr=yield_data["yieldRel_HMS_notrack"]*yield_data["uncern_yieldRel_HMS_notrack"],color='black',linestyle='None',zorder=3,label="_nolegend_")
        plt.scatter(yield_data["rate_HMS"]/1000,yield_data["yieldRel_HMS_notrack"],color='blue',zorder=4,label="_nolegend_")
        #yield_data["m0_rate_HMS_notrack"] = linear_plot(yield_data["rate_HMS"]/1000,yield_data["yieldRel_HMS_notrack"],None,yield_data["uncern_yieldRel_HMS_notrack"],xvalmax=max((yield_data["rate_HMS"])/1000)+5)
        #plt.errorbar(yield_data["rate_HMS"]/1000,yield_data["yieldRel_HMS_CPULT_notrack"],yerr=yield_data["yieldRel_HMS_CPULT_notrack"]*yield_data["uncern_yieldRel_HMS_CPULT_notrack"],color='black',linestyle='None',zorder=5)
        #plt.scatter(yield_data["rate_HMS"]/1000,yield_data["yieldRel_HMS_CPULT_notrack"],color='red',zorder=6)
        
        if(np.isnan(yield_data["yieldRel_HMS_notrack"][0])== False):
            yfit = (slope_HMS_ntrVSrate)*(yield_data["rate_HMS"]/1000) + inter_HMS_ntrVSrate
            plt.plot(yield_data["rate_HMS"]/1000, yfit, color = 'green', label = "slope = %f +/- %f" %(slope_HMS_ntrVSrate, d_slope_HMS_ntrVSrate) + "\n intercept = %f +/- %f" %(inter_HMS_ntrVSrate, d_inter_HMS_ntrVSrate))
    
        
        
        plt.ylabel('Rel. Yield no track', fontsize=RatePlotFontSize)
        plt.xlabel('HMS %s Rate [kHz]' % (str(HMSscaler)), fontsize=RatePlotFontSize)
        plt.legend(prop={'size' :7})
        if target == 'LD2' :
            plt.title('HMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize=RatePlotFontSize)
        elif target == 'LH2' :
            plt.title('HMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize=RatePlotFontSize)
        else :
            plt.title('HMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize=RatePlotFontSize)
    
        #HMS plot track
        plt.subplot(2,3,3)    
        plt.grid(zorder=1)
        plt.xlim(0,(max(yield_data["rate_HMS"])/1000)+5)
        plt.ylim(0.9,1.1)
        plt.plot([0,(max(yield_data["rate_HMS"])/1000)+5], [1,1], 'r-',zorder=2)
        
        if(np.isnan(yield_data["yieldRel_HMS_track"][0])== False):
            a_fit_HMS_trVSrate, cov_HMS_trVSrate = curve_fit(linfunc,yield_data["rate_HMS"]/1000, yield_data["yieldRel_HMS_track"], sigma=yield_data["uncern_yieldRel_HMS_track"], absolute_sigma = True)
            inter_HMS_trVSrate = a_fit_HMS_trVSrate[0]
            slope_HMS_trVSrate = a_fit_HMS_trVSrate[1]
            d_inter_HMS_trVSrate = np.sqrt(cov_HMS_trVSrate[0][0])
            d_slope_HMS_trVSrate = np.sqrt(cov_HMS_trVSrate[1][1])
        
        plt.errorbar(yield_data["rate_HMS"]/1000,yield_data["yieldRel_HMS_track"],yerr=yield_data["yieldRel_HMS_track"]*yield_data["uncern_yieldRel_HMS_track"],color='black',linestyle='None',zorder=3,label="_nolegend_")
        plt.scatter(yield_data["rate_HMS"]/1000,yield_data["yieldRel_HMS_track"],color='blue',zorder=4,label="_nolegend_")
       
        if(np.isnan(yield_data["yieldRel_HMS_track"][0])== False):
            yfit = (slope_HMS_trVSrate)*(yield_data["rate_HMS"]/1000) + inter_HMS_trVSrate
            plt.plot(yield_data["rate_HMS"]/1000, yfit, color = 'green', label = "slope = %f +/- %f" %(slope_HMS_trVSrate, d_slope_HMS_trVSrate) + "\n intercept = %f +/- %f" %(inter_HMS_trVSrate, d_inter_HMS_trVSrate))
    
        
        #yield_data["m0_rate_HMS_track"] = linear_plot(yield_data["rate_HMS"]/1000,yield_data["yieldRel_HMS_track"],None,yield_data["uncern_yieldRel_HMS_track"],xvalmax=max((yield_data["rate_HMS"])/1000)+5)
        #plt.errorbar(yield_data["rate_HMS"]/1000,yield_data["yieldRel_HMS_CPULT_track"],yerr=yield_data["yieldRel_HMS_CPULT_track"]*yield_data["uncern_yieldRel_HMS_CPULT_track"],color='black',linestyle='None',zorder=5)
        #plt.scatter(yield_data["rate_HMS"]/1000,yield_data["yieldRel_HMS_CPULT_track"],color='red',zorder=6)
        plt.ylabel('Rel. Yield track', fontsize=RatePlotFontSize)
        plt.xlabel('HMS %s Rate [kHz]' % (str(HMSscaler)), fontsize=RatePlotFontSize)
        plt.legend(prop={'size' :7})
        if target == 'LD2' :
            plt.title('HMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize=RatePlotFontSize)
        elif target == 'LH2' :
            plt.title('HMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize=RatePlotFontSize)
        else :
            plt.title('HMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize=RatePlotFontSize)

    if(not NO_SHMS):
        #SHMS plot scaler
        plt.subplot(3,3,4)    
        plt.grid(zorder=1)
        plt.xlim(0,(max(yield_data["rate_SHMS"])/1000)+5)
        plt.ylim(0.9,1.1)
        plt.plot([0,(max(yield_data["rate_SHMS"])/1000)+5], [1,1], 'r-',zorder=2)
        plt.errorbar(yield_data["rate_SHMS"]/1000,yield_data["yieldRel_SHMS_scaler"],yerr=yield_data["yieldRel_SHMS_scaler"]*yield_data["uncern_yieldRel_SHMS_scaler"],color='black',linestyle='None',zorder=3,label="_nolegend_")
        plt.scatter(yield_data["rate_SHMS"]/1000,yield_data["yieldRel_SHMS_scaler"],color='blue',zorder=4,label="_nolegend_")
        #yield_data["m0_rate_SHMS_scaler"] = linear_plot(yield_data["rate_SHMS"]/1000,yield_data["yieldRel_SHMS_scaler"],None,yield_data["uncern_yieldRel_SHMS_scaler"],xvalmax=max((yield_data["rate_SHMS"])/1000)+5)
    
        plt.ylabel('Rel. Yield %s' % (str(SHMSscaler)), fontsize=RatePlotFontSize)
        plt.xlabel('SHMS %s Rate [kHz]' % (str(SHMSscaler)), fontsize =RatePlotFontSize)
        plt.legend(prop={'size' :7})
        if target == 'LD2' :
            plt.title('SHMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =RatePlotFontSize)
        elif target == 'LH2' :
            plt.title('SHMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =RatePlotFontSize)
        else :
            plt.title('SHMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =RatePlotFontSize)
    
        #SHMS plot no track
        plt.subplot(3,3,5)    
        plt.grid(zorder=1)
        plt.xlim(0,(max(yield_data["rate_SHMS"])/1000)+5)
        plt.ylim(0.9,1.1)
        plt.plot([0,(max(yield_data["rate_SHMS"])/1000)+5], [1,1], 'r-',zorder=2)
        
        plt.errorbar(yield_data["rate_SHMS"]/1000,yield_data["yieldRel_SHMS_notrack"],yerr=yield_data["yieldRel_SHMS_notrack"]*yield_data["uncern_yieldRel_SHMS_notrack"],color='black',linestyle='None',zorder=3,label="_nolegend_")
        plt.scatter(yield_data["rate_SHMS"]/1000,yield_data["yieldRel_SHMS_notrack"],color='blue',zorder=4,label="_nolegend_")
        #yield_data["m0_rate_SHMS_notrack"] = linear_plot(yield_data["rate_SHMS"]/1000,yield_data["yieldRel_SHMS_notrack"],None,yield_data["uncern_yieldRel_SHMS_notrack"],xvalmax=max((yield_data["rate_SHMS"])/1000)+5)
        #plt.errorbar(yield_data["rate_SHMS"]/1000,yield_data["yieldRel_SHMS_CPULT_notrack"],yerr=yield_data["yieldRel_SHMS_CPULT_notrack"]*yield_data["uncern_yieldRel_SHMS_CPULT_notrack"],color='black',linestyle='None',zorder=5)
        #plt.scatter(yield_data["rate_SHMS"]/1000,yield_data["yieldRel_SHMS_CPULT_notrack"],color='red',zorder=6)
        plt.ylabel('Rel. Yield no track', fontsize=RatePlotFontSize)
        plt.xlabel('SHMS %s Rate [kHz]' % (str(SHMSscaler)), fontsize =RatePlotFontSize)
        plt.legend(prop={'size' :7})
        if target == 'LD2' :
            plt.title('SHMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =RatePlotFontSize)
        elif target == 'LH2' :
            plt.title('SHMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =RatePlotFontSize)
        else :
            plt.title('SHMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =RatePlotFontSize)

        #SHMS plot track
        plt.subplot(3,3,6)    
        plt.grid(zorder=1)
        plt.xlim(0,(max(yield_data["rate_SHMS"])/1000)+5)
        plt.ylim(0.9,1.1)
        plt.plot([0,(max(yield_data["rate_SHMS"])/1000)+5], [1,1], 'r-',zorder=2)
        
        if(np.isnan(yield_data["yieldRel_SHMS_track"][0])== False):
            a_fit_SHMS_trVSrate, cov_SHMS_trVSrate = curve_fit(linfunc,yield_data["rate_SHMS"]/1000, yield_data["yieldRel_SHMS_track"], sigma=yield_data["uncern_yieldRel_SHMS_track"], absolute_sigma = True)
            inter_SHMS_trVSrate = a_fit_SHMS_trVSrate[0]
            slope_SHMS_trVSrate = a_fit_SHMS_trVSrate[1]
            d_inter_SHMS_trVSrate = np.sqrt(cov_SHMS_trVSrate[0][0])
            d_slope_SHMS_trVSrate = np.sqrt(cov_SHMS_trVSrate[1][1])
        
        plt.errorbar(yield_data["rate_SHMS"]/1000,yield_data["yieldRel_SHMS_track"],yerr=yield_data["yieldRel_SHMS_track"]*yield_data["uncern_yieldRel_SHMS_track"],color='black',linestyle='None',zorder=3,label="_nolegend_")
        plt.scatter(yield_data["rate_SHMS"]/1000,yield_data["yieldRel_SHMS_track"],color='blue',zorder=4,label="_nolegend_")
        #yield_data["m0_rate_SHMS_track"] = linear_plot(yield_data["rate_SHMS"]/1000,yield_data["yieldRel_SHMS_track"],None,yield_data["uncern_yieldRel_SHMS_track"],xvalmax=max((yield_data["rate_SHMS"])/1000)+5)
        #plt.errorbar(yield_data["rate_SHMS"]/1000,yield_data["yieldRel_SHMS_CPULT_track"],yerr=yield_data["yieldRel_SHMS_CPULT_track"]*yield_data["uncern_yieldRel_SHMS_CPULT_track"],color='black',linestyle='None',zorder=5)
        #plt.scatter(yield_data["rate_SHMS"]/1000,yield_data["yieldRel_SHMS_CPULT_track"],color='red',zorder=6)
        
        if(np.isnan(yield_data["yieldRel_SHMS_track"][0])== False):
            yfit = (slope_SHMS_trVSrate)*(yield_data["rate_SHMS"]/1000) + inter_SHMS_trVSrate
            plt.plot(yield_data["rate_SHMS"]/1000, yfit, color = 'green', label = "slope = %f +/- %f" %(slope_SHMS_trVSrate, d_slope_SHMS_trVSrate) + "\n intercept = %f +/- %f" %(inter_SHMS_trVSrate, d_inter_SHMS_trVSrate))
    
        
        plt.ylabel('Rel. Yield track', fontsize=RatePlotFontSize)
        plt.xlabel('SHMS %s Rate [kHz]' % (str(SHMSscaler)), fontsize =RatePlotFontSize)
        plt.legend(prop={'size' :7})
        if target == 'LD2' :
            plt.title('SHMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =RatePlotFontSize)
        elif target == 'LH2' :
            plt.title('SHMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =RatePlotFontSize)
        else :
            plt.title('SHMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =RatePlotFontSize)

    ## Begin Coin ##
    ## Scaler ##
    plt.subplot(3,3,7)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.ylim(0.9,3)
    plt.plot([0,100], [1,1], 'r-',zorder=2)
    plt.errorbar(yield_data["rate_COIN"]/1000,yield_data["yieldRel_COIN_scaler"],yerr=yield_data["yieldRel_COIN_scaler"]*yield_data["uncern_yieldRel_COIN_scaler"],color='black',linestyle='None',zorder=3,label="_nolegend_")
    plt.scatter(yield_data["rate_COIN"]/1000,yield_data["yieldRel_COIN_scaler"],color='blue',zorder=4,label="_nolegend_")
    
    plt.ylabel('Rel. Yield %s' % (str(COINscaler)), fontsize=RatePlotFontSize)
    plt.xlabel('COIN %s [kHz]' % (str(COINscaler)), fontsize =RatePlotFontSize)
    plt.legend(prop={'size' :7})    
    if target == 'LD2' :
        plt.title('COIN %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =RatePlotFontSize)
    elif target == 'LH2' :
        plt.title('COIN %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =RatePlotFontSize)
    else :
        plt.title('COIN %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =RatePlotFontSize)
            
    ## No Track ##        
    plt.subplot(3,3,8)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.ylim(0.9,1.2)
    plt.plot([0,100], [1,1], 'r-',zorder=2)
    plt.errorbar(yield_data["rate_COIN"]/1000,yield_data["yieldRel_COIN_notrack"],yerr=yield_data["yieldRel_COIN_notrack"]*yield_data["uncern_yieldRel_COIN_notrack"],color='black',linestyle='None',zorder=3,label="_nolegend_")
    plt.scatter(yield_data["rate_COIN"]/1000,yield_data["yieldRel_COIN_notrack"],color='blue',zorder=4,label="_nolegend_")
    
    plt.ylabel('Rel. Yield no track', fontsize=RatePlotFontSize)
    plt.xlabel('COIN %s [kHz]' % (str(COINscaler)), fontsize =RatePlotFontSize)
    plt.legend(prop={'size' :7})    
    if target == 'LD2' :
        plt.title('COIN %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =RatePlotFontSize)
    elif target == 'LH2' :
        plt.title('COIN %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =RatePlotFontSize)
    else :
        plt.title('COIN %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =RatePlotFontSize)
    
    ## Track ##
    plt.subplot(3,3,9)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.ylim(0.9,1.2)
    plt.plot([0,100], [1,1], 'r-',zorder=2)
    plt.errorbar(yield_data["rate_COIN"]/1000,yield_data["yieldRel_COIN_track"],yerr=yield_data["yieldRel_COIN_track"]*yield_data["uncern_yieldRel_COIN_track"],color='black',linestyle='None',zorder=3,label="_nolegend_")
    plt.scatter(yield_data["rate_COIN"]/1000,yield_data["yieldRel_COIN_track"],color='blue',zorder=4,label="_nolegend_")
    
    plt.ylabel('Rel. Yield track', fontsize=RatePlotFontSize)
    plt.xlabel('COIN %s [kHz]' % (str(COINscaler)), fontsize =RatePlotFontSize)
    plt.legend(prop={'size' :7})    
    if target == 'LD2' :
        plt.title('COIN %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =RatePlotFontSize)
    elif target == 'LH2' :
        plt.title('COIN %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =RatePlotFontSize)
    else :
        plt.title('COIN %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =RatePlotFontSize)
    

    #plt.tight_layout()
    plt.savefig(SCRIPTPATH+'/luminosity/OUTPUTS/plots/Yield_%s_%s.png' % (out_f.split("yield_data_")[1].replace(".csv",""),"raterelYieldPlot"))
            

    #########################################################################################################################################################

    edtmPlot = plt.figure(figsize=(12,8))

    #Ratio accp/total scaler vs current
    plt.subplot(2,3,1)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.scatter(yield_data["current"],(yield_data["HMS_scaler_accp"]+yield_data["SHMS_scaler_accp"])/(yield_data["HMSTRIG_scaler"]+yield_data["SHMSTRIG_scaler"]),color='blue',zorder=4,label="_nolegend_")
    plt.ylabel('Accp/Total scaler count', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =16)
    if target == 'LD2' :
        plt.title('LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)

    #Scaler EDTM rate vs current
    plt.subplot(2,3,2)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.scatter(yield_data["current"],yield_data["sent_edtm"]/yield_data["time"],color='blue',zorder=4,label="_nolegend_")
    plt.scatter(yield_data["current"],yield_data["sent_edtm_PS"]/yield_data["time"],color='red',zorder=4,label="_nolegend_")
    plt.ylabel('Scaler EDTM rate [Hz]', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =16)
    if target == 'LD2' :
        plt.title('HMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('HMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('HMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)

    #EDTM vs HMS Rate
    plt.subplot(2,3,3)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.scatter(yield_data["rate_HMS"]/1000,yield_data["accp_edtm"]/(yield_data["time"]*1000),color='blue',zorder=4,label="_nolegend_")
    plt.ylabel('EDTM Rate [kHz]', fontsize=16)
    plt.xlabel('HMS %s Rate [kHz]' % (str(HMSscaler)), fontsize =16)
    if target == 'LD2' :
        plt.title('HMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('HMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('HMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)

    #TLT vs Current
    plt.subplot(2,3,4)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.errorbar(yield_data["current"],yield_data["TLT"],yerr=yield_data["TLT"]*yield_data["uncern_TLT"],color='black',linestyle='None',zorder=3,label="_nolegend_")
    plt.scatter(yield_data["current"],yield_data["TLT"],color='blue',zorder=4,label="_nolegend_")
    plt.ylabel('TLT', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =16)
    if target == 'LD2' :
        plt.title('LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)

    #Time vs current
    plt.subplot(2,3,5)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.scatter(yield_data["current"],yield_data["time"]/70,color='blue',zorder=4,label="_nolegend_")
    plt.ylabel('Time [min]', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =16)
    if target == 'LD2' :
        plt.title('HMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('HMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('HMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)


    #EDTM vs SHMS Rate
    plt.subplot(2,3,6)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.scatter(yield_data["rate_SHMS"]/1000,yield_data["accp_edtm"]/(yield_data["time"]*1000),color='blue',zorder=4,label="_nolegend_")
    plt.ylabel('EDTM Rate [kHz]', fontsize=16)
    plt.xlabel('SHMS %s Rate [kHz]' % (str(SHMSscaler)), fontsize =16)
    if target == 'LD2' :
        plt.title('SHMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('SHMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('SHMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)

    plt.tight_layout()             
    plt.savefig(SCRIPTPATH+'/luminosity/OUTPUTS/plots/Yield_%s_%s.png' % (out_f.split("yield_data_")[1].replace(".csv",""),"edtmPlot"))
            
    #########################################################################################################################################################

    logPlot = plt.figure(figsize=(12,8))

    if(not NO_HMS):
        #HMS plot scaler
        plt.subplot(2,4,1)    
        plt.grid(zorder=1)
        plt.xlim(0,100)
        plt.ylim(0.9,1.04)
        plt.plot([0,100], [1,1], 'r-',zorder=2)
        
        ''' if(np.isnan(yield_data["yieldRel_HMS_scaler"][0])== False):
            a_fit_HMS_scalerVScurrent, cov_HMS_scalerVScurrent = curve_fit(linfunc,yield_data["current"], yield_data["yieldRel_HMS_scaler"], sigma=yield_data["uncern_yieldRel_HMS_scaler"], absolute_sigma = True)
            inter_HMS_scalerVScurrent = a_fit_HMS_scalerVScurrent[0]
            slope_HMS_scalerVScurrent = a_fit_HMS_scalerVScurrent[1]
            d_inter_HMS_scalerVScurrent = np.sqrt(cov_HMS_scalerVScurrent[0][0])
            d_slope_HMS_scalerVScurrent = np.sqrt(cov_HMS_scalerVScurrent[1][1])'''
        
        plt.errorbar(yield_data["current"],yield_data["yieldRel_HMS_scaler"],yerr=yield_data["yieldRel_HMS_scaler"]*yield_data["uncern_yieldRel_HMS_scaler"],color='black',linestyle='None',zorder=3,label="_nolegend_")
        plt.scatter(yield_data["current"],yield_data["yieldRel_HMS_scaler"],color='blue',zorder=4, label = "_nolegend_")
        
        '''if(np.isnan(yield_data["yieldRel_HMS_scaler"][0])== False):
            yfit = (slope_HMS_scalerVScurrent/inter_HMS_scalerVScurrent)*yield_data["current"] + 1.0
            plt.plot(yield_data["current"], yfit, color = 'green', label = "slope = %f +/- %f" %(slope_HMS_scalerVScurrent, d_slope_HMS_scalerVScurrent) + "\n intercept = %f +/- %f" %(inter_HMS_scalerVScurrent, d_inter_HMS_scalerVScurrent))
        '''
        
        plt.ylabel('Rel. Yield %s' % (str(HMSscaler)), fontsize=16)
        plt.xlabel('Current [uA]', fontsize =12)
        plt.legend(prop={'size' :7})
        if target == 'LD2' :
            plt.title('HMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
        elif target == 'LH2' :
            plt.title('HMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
        else :
            plt.title('HMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    
        #HMS plot track
        plt.subplot(2,4,2)    
        plt.grid(zorder=1)
        #plt.xlim(0,100)
        plt.ylim(0.9,1.04)
        plt.plot([0,100], [1,1], 'r-',zorder=2)
        
        '''if(np.isnan(yield_data["yieldRel_HMS_track"][0])== False):
            a_fit_HMS_trVScurrent, cov_HMS_trVScurrent = curve_fit(linfunc,yield_data["current"], yield_data["yieldRel_HMS_track"], sigma=yield_data["uncern_yieldRel_HMS_track"], absolute_sigma = True)
            inter_HMS_trVScurrent = a_fit_HMS_trVScurrent[0]
            slope_HMS_trVScurrent = a_fit_HMS_trVScurrent[1]
            d_inter_HMS_trVScurrent = np.sqrt(cov_HMS_trVScurrent[0][0])
            d_slope_HMS_trVScurrent = np.sqrt(cov_HMS_trVScurrent[1][1])'''
        
        plt.errorbar(yield_data["current"],yield_data["yieldRel_HMS_track"],yerr=yield_data["yieldRel_HMS_track"]*yield_data["uncern_yieldRel_HMS_track"],color='black',linestyle='None',zorder=3,label="_nolegend_")
        plt.scatter(yield_data["current"],yield_data["yieldRel_HMS_track"],color='blue',zorder=4,label="_nolegend_")
        
        '''if(np.isnan(yield_data["yieldRel_HMS_track"][0])== False):
            yfit = (slope_HMS_trVScurrent/inter_HMS_trVScurrent)*yield_data["current"] + 1.0
            plt.plot(yield_data["current"], yfit, color = 'green', label = "slope = %f +/- %f" %(slope_HMS_trVScurrent, d_slope_HMS_trVScurrent) + "\n intercept = %f +/- %f" %(inter_HMS_trVScurrent, d_inter_HMS_trVScurrent))
         '''
        
        plt.ylabel('Rel. Yield track', fontsize=16)
        plt.xlabel('Current [uA]', fontsize =12)
        plt.legend(prop={'size' :7})
        if target == 'LD2' :
            plt.title('HMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
        elif target == 'LH2' :
            plt.title('HMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
        else :
            plt.title('HMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    
        #HMS track eff vs HMS Rate
        plt.subplot(2,4,3)    
        plt.grid(zorder=1)
        #plt.xlim(0,100)
        plt.errorbar(yield_data["rate_HMS"]/1000,yield_data["HMS_track"],yerr=yield_data["HMS_track"]*yield_data["HMS_track_uncern"],color='black',linestyle='None',zorder=3,label="_nolegend_")
        plt.scatter(yield_data["rate_HMS"]/1000,yield_data["HMS_track"],color='blue',zorder=4,label="_nolegend_")
        plt.ylabel('HMS track eff', fontsize=16)
        plt.xlabel('HMS %s Rate [kHz]' % (str(HMSscaler)), fontsize =12)
        if target == 'LD2' :
            plt.title('HMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
        elif target == 'LH2' :
            plt.title('HMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
        else :
            plt.title('HMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)

    #TLT vs Current 
    '''plt.subplot(2,4,4)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.errorbar(yield_data["current"],yield_data["TLT"],yerr=yield_data["TLT"]*yield_data["uncern_TLT"],color='black',linestyle='None',zorder=3,label="_nolegend_")
    plt.scatter(yield_data["current"],yield_data["TLT"],color='blue',zorder=4,label="_nolegend_")
    plt.ylabel('TLT', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =12)
    if target == 'LD2' :
        plt.title('LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    elif target == 'LH2' :
        plt.title('LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    else :
        plt.title('Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    '''

    #TLT vs Rate
    plt.subplot(2,4,4)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.ylabel('TLT', fontsize=16)
    
    plt.xlabel('Total Trigger Rate (kHz)', fontsize =12)
    plt.errorbar(yield_data["rate_TOTAL"]/1000,yield_data["TLT"],yerr=yield_data["TLT"]*yield_data["uncern_TLT"],color='black',linestyle='None',zorder=3,label="_nolegend_")
    plt.scatter(yield_data["rate_TOTAL"]/1000,yield_data["TLT"],color='red',zorder=4,label="_nolegend_") #TLT using EDTM
    
    plt.errorbar(yield_data["rate_TOTAL"]/1000,yield_data["TLT_ELT"],yerr=yield_data["TLT_ELT"]*np.sqrt(yield_data["uncern_TLT_ELT"]),color='black',linestyle='None',zorder=3,label="_nolegend_")
    plt.scatter(yield_data["rate_TOTAL"]/1000,yield_data["TLT_ELT"],color='blue',zorder=4,label="_nolegend_") #TLT using CPULT
       
    plt.legend(prop={'size' :7})
    
    
    
    if target == 'LD2' :
        plt.title('LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    elif target == 'LH2' :
        plt.title('LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    else :
        plt.title('Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)



    #CPULT vs Current
    plt.subplot(2,4,8)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.errorbar(yield_data["current"],yield_data["CPULT_phys"],yerr=yield_data["CPULT_phys"]*yield_data["uncern_CPULT_phys"],color='black',linestyle='None',zorder=3,label="_nolegend_")
    plt.scatter(yield_data["current"],yield_data["CPULT_phys"],color='blue',zorder=4,label="_nolegend_")
    plt.ylabel('CPULT', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =12)
    if target == 'LD2' :
        plt.title('LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    elif target == 'LH2' :
        plt.title('LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    else :
        plt.title('Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)

    if(not NO_SHMS):
        #SHMS plot scaler
        plt.subplot(2,4,5)    
        plt.grid(zorder=1)
        plt.xlim(0,100)
        plt.ylim(0.9,1.04)
        plt.plot([0,100], [1,1], 'r-',zorder=2)
        
        '''if(np.isnan(yield_data["yieldRel_SHMS_scaler"][0])== False):
            a_fit_SHMS_scalerVScurrent, cov_SHMS_scalerVScurrent = curve_fit(linfunc,yield_data["current"], yield_data["yieldRel_SHMS_scaler"], sigma=yield_data["uncern_yieldRel_SHMS_scaler"], absolute_sigma = True)
            inter_SHMS_scalerVScurrent = a_fit_SHMS_scalerVScurrent[0]
            slope_SHMS_scalerVScurrent = a_fit_SHMS_scalerVScurrent[1]
            d_inter_SHMS_scalerVScurrent = np.sqrt(cov_SHMS_scalerVScurrent[0][0])
            d_slope_SHMS_scalerVScurrent = np.sqrt(cov_SHMS_scalerVScurrent[1][1])'''
        
        plt.errorbar(yield_data["current"],yield_data["yieldRel_SHMS_scaler"],yerr=yield_data["yieldRel_SHMS_scaler"]*yield_data["uncern_yieldRel_SHMS_scaler"],color='black',linestyle='None',zorder=3,label="_nolegend_")
        plt.scatter(yield_data["current"],yield_data["yieldRel_SHMS_scaler"],color='blue',zorder=4,label="_nolegend_")
        
        '''if(np.isnan(yield_data["yieldRel_SHMS_scaler"][0])== False):
            yfit = (slope_SHMS_scalerVScurrent/inter_SHMS_scalerVScurrent)*yield_data["current"] + 1.0
            plt.plot(yield_data["current"], yfit, color = 'green', label = "slope = %f +/- %f" %(slope_SHMS_scalerVScurrent, d_slope_SHMS_scalerVScurrent) + "\n intercept = %f +/- %f" %(inter_SHMS_scalerVScurrent, d_inter_SHMS_scalerVScurrent))
        '''
        
        plt.ylabel('Rel. Yield %s' % (str(SHMSscaler)), fontsize=16)
        plt.xlabel('Current [uA]', fontsize =12)
        plt.legend(prop={'size' :7})
        if target == 'LD2' :
            plt.title('SHMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
        elif target == 'LH2' :
            plt.title('SHMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
        else :
            plt.title('SHMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    
        #SHMS plot track
        plt.subplot(2,4,6)    
        plt.grid(zorder=1)
        #plt.xlim(0,100)
        plt.ylim(0.9,1.04)
        plt.plot([0,100], [1,1], 'r-',zorder=2)
        
        ''' if(np.isnan(yield_data["yieldRel_SHMS_track"][0])== False):
            a_fit_SHMS_trVScurrent, cov_SHMS_trVScurrent = curve_fit(linfunc,yield_data["current"], yield_data["yieldRel_SHMS_track"], sigma=yield_data["uncern_yieldRel_SHMS_track"], absolute_sigma = True)
            inter_SHMS_trVScurrent = a_fit_SHMS_trVScurrent[0]
            slope_SHMS_trVScurrent = a_fit_SHMS_trVScurrent[1]
            d_inter_SHMS_trVScurrent = np.sqrt(cov_SHMS_trVScurrent[0][0])
            d_slope_SHMS_trVScurrent = np.sqrt(cov_SHMS_trVScurrent[1][1])'''
        
        plt.errorbar(yield_data["current"],yield_data["yieldRel_SHMS_track"],yerr=yield_data["yieldRel_SHMS_track"]*yield_data["uncern_yieldRel_SHMS_track"],color='black',linestyle='None',zorder=3,label="_nolegend_")
        plt.scatter(yield_data["current"],yield_data["yieldRel_SHMS_track"],color='blue',zorder=4,label="_nolegend_")
        
        ''' if(np.isnan(yield_data["yieldRel_SHMS_track"][0])== False):
            yfit = (slope_SHMS_trVScurrent/inter_SHMS_trVScurrent)*yield_data["current"] + 1.0
            plt.plot(yield_data["current"], yfit, color = 'green', label = "slope = %f +/- %f" %(slope_SHMS_trVScurrent, d_slope_SHMS_trVScurrent) + "\n intercept = %f +/- %f" %(inter_SHMS_trVScurrent, d_inter_SHMS_trVScurrent))
        '''     
        plt.ylabel('Rel. Yield track', fontsize=16)
        plt.xlabel('Current [uA]', fontsize =12)
        plt.legend(prop={'size' :7})
        if target == 'LD2' :
            plt.title('SHMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
        elif target == 'LH2' :
            plt.title('SHMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
        else :
            plt.title('SHMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    
        #SHMS track vs SHMS Rate
        plt.subplot(2,4,7)    
        plt.grid(zorder=1)
        #plt.xlim(0,100)
        plt.errorbar(yield_data["rate_SHMS"]/1000,yield_data["SHMS_track"],yerr=yield_data["SHMS_track"]*yield_data["SHMS_track_uncern"],color='black',linestyle='None',zorder=3,label="_nolegend_")
        plt.scatter(yield_data["rate_SHMS"]/1000,yield_data["SHMS_track"],color='blue',zorder=4,label="_nolegend_")
        plt.ylabel('SHMS track eff', fontsize=16)
        plt.xlabel('SHMS %s Rate [kHz]' % (str(SHMSscaler)), fontsize =12)
        if target == 'LD2' :
            plt.title('SHMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
        elif target == 'LH2' :
            plt.title('SHMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
        else :
            plt.title('SHMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    
    plt.tight_layout()
    plt.savefig(SCRIPTPATH+'/luminosity/OUTPUTS/plots/Yield_%s_%s.png' % (out_f.split("yield_data_")[1].replace(".csv",""),"logPlot"))
            
    #plt.show()

    print("\nYield info...\n",yield_data[["run number","yieldRel_HMS_scaler","uncern_yieldRel_HMS_scaler","yieldRel_SHMS_scaler","uncern_yieldRel_SHMS_scaler","yieldRel_HMS_notrack","uncern_yieldRel_HMS_notrack","yieldRel_SHMS_notrack","uncern_yieldRel_SHMS_notrack","yieldRel_HMS_track","uncern_yieldRel_HMS_track","yieldRel_SHMS_track","uncern_yieldRel_SHMS_track"]])

    # x axis = current (uA)
    #print("HMS scaler v. current: slope + intercept = ", poly_HMS_scalerVScurrent,"\n" )
    #print("HMS notrack v. current: slope + intercept = " , poly_HMS_ntrVScurrent, "\n")
    #print("HMS track v. current: slope + intercept = ", poly_HMS_trVScurrent, "\n")
    
    #print("SHMS scaler v. current: slope + intercept = ", poly_SHMS_scalerVScurrent,"\n" )
    #print("SHMS notrack v. current: slope + intercept = " , poly_SHMS_ntrVScurrent, "\n")
    #print("SHMS track v. current: slope + intercept = ", poly_SHMS_trVScurrent, "\n")
    
    # x axis = rate (kHz)
    #print("HMS scaler v. rate: slope + intercept = ", poly_HMS_scalerVSrate, "\n")
    #print("HMS notrack v. rate: slope + intercept = ", poly_HMS_ntrVSrate, "\n")
    #print("HMS track v. rate: slope + intercept = ", poly_HMS_trVSrate, "\n")
    
    #print("SHMS scaler v. rate: slope + intercept = ", poly_SHMS_scalerVSrate, "\n")
    #print("SHMS notrack v. rate: slope + intercept = ", poly_SHMS_ntrVSrate, "\n")
    #print("SHMS track v. rate: slope + intercept = ", poly_SHMS_trVSrate, "\n")
    #print(slope_HMS_scalerVScurrent, d_slope_HMS_scalerVScurrent, inter_HMS_scalerVScurrent, d_inter_HMS_scalerVScurrent)
    
   
   
    return yield_data
    
################################################################################################################################################

def debug():
    '''
    Prints various values to debug, customize to your heart's content
    '''
    data = mergeDicts()
    print("\n\n=======================")
    print("DEBUG data")
    print("=======================")
    ### Debug prints
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        print(data[["run number","sent_edtm","sent_edtm_PS","sent_edtm_HMS","sent_edtm_SHMS","accp_edtm","TLT","CPULT_phys","charge","current","time","HMS_track","SHMS_track"]])
   # print("EDTM scaler rate: ", data["sent_edtm"]/data["time"])
   # print("Accepted EDTM rate: ", data["accp_edtm"]/data["time"])
   # print("Run numbers: ", data["run number"].sort_values())
   # print("HMS scaler ratio",data["HMS_scaler_accp"]/data["HMSTRIG_scaler"])
   # print("SHMS scaler ratio",data["SHMS_scaler_accp"]/data["SHMSTRIG_scaler"])

    print("HMS Cal etotnorm\n",data[["h_int_etotnorm_evts","current"]])
    print("SHMS Cal etotnorm\n",data[["p_int_etotnorm_evts","current"]])
    print("HMS yield no track\n",data["yield_HMS_notrack"])
    print("SHMS yield no track\n",data["yield_SHMS_notrack"])
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        print("\n\nHMS track uncertanties\n",data[["current","uncern_HMS_evts_notrack", "uncern_TLT", "uncern_CPULT_phys", "HMS_track_uncern"]])
        print("\n\nSHMS track uncertanties\n",data[["current","uncern_SHMS_evts_notrack", "uncern_TLT", "uncern_CPULT_phys", "SHMS_track_uncern"]])
    print("\n\nHMS yield track uncern\n",data["uncern_yieldRel_HMS_track"])
    print("SHMS yield track uncern\n",data["uncern_yieldRel_SHMS_track"])
    ###
    print("=======================\n\n")

    print("Average HMS",np.average(data["yieldRel_HMS_scaler"]))
    print("Average SHMS",np.average(data["yieldRel_SHMS_scaler"]))
    print("Average I",np.average(data["current"]))
    
   
    
    

################################################################################################################################################

def main():


    debug()

    # Plot yields
    yield_data = plot_yield()

    table = mergeDicts()
    
    file_exists = os.path.isfile(out_f)

    table.to_csv(out_f, index=False, header=True, mode='w+',)
    #print(yield_data["run number", "yieldRel_COIN_track", "uncern_yieldRel_COIN_track"],file=out_f) 



if __name__ == '__main__':
    main()
