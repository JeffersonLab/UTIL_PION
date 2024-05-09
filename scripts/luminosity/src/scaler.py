#! /usr/bin/python
#
# Description: This is where the scaler variables for the yield calculations are formulated.
# Variables calculated: SHMS_PS, HMS_PS, time, charge, SHMSTRIG_scaler, HMSTRIG_scaler, CPULT_scaler, CPULT_scaler_uncern, HMS_eLT, HMS_eLT_uncern, SHMS_eLT, SHMS_eLT_uncern, sent_edtm
# ================================================================
# Time-stamp: "2023-06-28 10:51:41 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import uproot as up
import numpy as np
import math

################################################################################################################################################

def scaler(PS_names, HMS_PS, SHMS_PS, COIN_PS, thres_curr, report_current, runNum, MaxEvent, s_tree):

    '''
    SCALER TREE, TSP
    '''

    # s_evts = len(s_tree)
    if SHMS_PS == None:
        s_evts = s_tree.array("H.BCM4A.scaler")

        P_BCM4A_scalerCharge = s_tree.array("H.BCM4A.scalerCharge")
        P_BCM2_scalerCharge = s_tree.array("H.BCM2.scalerCharge")
        P_BCM4B_scalerCharge = s_tree.array("H.BCM4B.scalerCharge")
        P_BCM1_scalerCharge = s_tree.array("H.BCM1.scalerCharge")
        P_BCM4C_scalerCharge = s_tree.array("H.BCM4C.scalerCharge")

        P_BCM4A_scalerCurrent = s_tree.array("H.BCM4A.scalerCurrent")
        P_BCM2_scalerCurrent = s_tree.array("H.BCM2.scalerCurrent")
        P_BCM4B_scalerCurrent = s_tree.array("H.BCM4B.scalerCurrent")
        P_BCM1_scalerCurrent = s_tree.array("H.BCM1.scalerCurrent")
        P_BCM4C_scalerCurrent = s_tree.array("H.BCM4C.scalerCurrent")

        P_1Mhz_scalerTime = s_tree.array("H.1MHz.scalerTime")

        P_pTRIG1_scaler = s_tree.array("H.pTRIG1.scaler") #heinricn - changed from hTRIG to pTRIG because hTRIG isn't filled for some reason!
        P_pTRIG2_scaler = s_tree.array("H.pTRIG2.scaler")
        P_pTRIG3_scaler = s_tree.array("H.pTRIG3.scaler")
        P_pTRIG4_scaler = s_tree.array("H.pTRIG4.scaler")
        P_pTRIG5_scaler = s_tree.array("H.pTRIG5.scaler")
        P_pTRIG6_scaler = s_tree.array("H.pTRIG6.scaler")

        P_pL1ACCP_scaler = s_tree.array("H.hL1ACCP.scaler")
        P_pPRE40_scaler = s_tree.array("H.hPRE40.scaler")
        P_pPRE100_scaler = s_tree.array("H.hPRE100.scaler")
        P_pPRE150_scaler = s_tree.array("H.hPRE150.scaler")
        P_pPRE200_scaler = s_tree.array("H.hPRE200.scaler")
        P_pPRE40_scaler = s_tree.array("H.hPRE40.scaler")
        P_pPRE100_scaler = s_tree.array("H.hPRE100.scaler")
        P_pPRE150_scaler = s_tree.array("H.hPRE150.scaler")
        P_pPRE200_scaler = s_tree.array("H.hPRE200.scaler")

        P_pEL_LO_LO_scaler = s_tree.array("H.hEL_LO_LO.scaler")
        P_pEL_LO_scaler = s_tree.array("H.hEL_LO.scaler")
        P_pEL_HI_scaler = s_tree.array("H.hEL_HI.scaler")
        P_pEL_REAL_scaler = s_tree.array("H.hEL_REAL.scaler")
        P_pEL_CLEAN_scaler = s_tree.array("H.hEL_CLEAN.scaler")
        P_pSTOF_scaler = s_tree.array("H.hSTOF.scaler")

        P_pEL_LO_LO_scaler = s_tree.array("H.hEL_LO_LO.scaler")
        P_pEL_LO_scaler = s_tree.array("H.hEL_LO.scaler")
        P_pEL_HI_scaler = s_tree.array("H.hEL_HI.scaler")
        P_pEL_REAL_scaler = s_tree.array("H.hEL_REAL.scaler")
        P_pEL_CLEAN_scaler = s_tree.array("H.hEL_CLEAN.scaler")
        P_pSTOF_scaler = s_tree.array("H.hSTOF.scaler")
        P_pPRHI_scaler = s_tree.array("H.PRHI.scaler")
        P_pPRLO_scaler = s_tree.array("H.PRLO.scaler")

        P_EDTM_scaler = s_tree.array("H.EDTM.scaler")

    else:
        s_evts = s_tree.array("P.BCM4A.scaler")

        P_BCM4A_scalerCharge = s_tree.array("P.BCM4A.scalerCharge")
        P_BCM2_scalerCharge = s_tree.array("P.BCM2.scalerCharge")
        P_BCM4B_scalerCharge = s_tree.array("P.BCM4B.scalerCharge")
        P_BCM1_scalerCharge = s_tree.array("P.BCM1.scalerCharge")
        P_BCM4C_scalerCharge = s_tree.array("P.BCM4C.scalerCharge")

        P_BCM4A_scalerCurrent = s_tree.array("P.BCM4A.scalerCurrent")
        P_BCM2_scalerCurrent = s_tree.array("P.BCM2.scalerCurrent")
        P_BCM4B_scalerCurrent = s_tree.array("P.BCM4B.scalerCurrent")
        P_BCM1_scalerCurrent = s_tree.array("P.BCM1.scalerCurrent")
        P_BCM4C_scalerCurrent = s_tree.array("P.BCM4C.scalerCurrent")

        P_1Mhz_scalerTime = s_tree.array("P.1MHz.scalerTime")

        P_pTRIG1_scaler = s_tree.array("P.pTRIG1.scaler")
        P_pTRIG2_scaler = s_tree.array("P.pTRIG2.scaler")
        P_pTRIG3_scaler = s_tree.array("P.pTRIG3.scaler")
        P_pTRIG4_scaler = s_tree.array("P.pTRIG4.scaler")
        P_pTRIG5_scaler = s_tree.array("P.pTRIG5.scaler")
        P_pTRIG6_scaler = s_tree.array("P.pTRIG6.scaler")

        P_pL1ACCP_scaler = s_tree.array("P.pL1ACCP.scaler")
        P_pPRE40_scaler = s_tree.array("P.pPRE40.scaler")
        P_pPRE100_scaler = s_tree.array("P.pPRE100.scaler")
        P_pPRE150_scaler = s_tree.array("P.pPRE150.scaler")
        P_pPRE200_scaler = s_tree.array("P.pPRE200.scaler")
        P_pPRE40_scaler = s_tree.array("P.pPRE40.scaler")
        P_pPRE100_scaler = s_tree.array("P.pPRE100.scaler")
        P_pPRE150_scaler = s_tree.array("P.pPRE150.scaler")
        P_pPRE200_scaler = s_tree.array("P.pPRE200.scaler")

        P_pEL_LO_LO_scaler = s_tree.array("P.pEL_LO_LO.scaler")
        P_pEL_LO_scaler = s_tree.array("P.pEL_LO.scaler")
        P_pEL_HI_scaler = s_tree.array("P.pEL_HI.scaler")
        P_pEL_REAL_scaler = s_tree.array("P.pEL_REAL.scaler")
        P_pEL_CLEAN_scaler = s_tree.array("P.pEL_CLEAN.scaler")
        P_pSTOF_scaler = s_tree.array("P.pSTOF.scaler")

        P_pEL_LO_LO_scaler = s_tree.array("P.pEL_LO_LO.scaler")
        P_pEL_LO_scaler = s_tree.array("P.pEL_LO.scaler")
        P_pEL_HI_scaler = s_tree.array("P.pEL_HI.scaler")
        P_pEL_REAL_scaler = s_tree.array("P.pEL_REAL.scaler")
        P_pEL_CLEAN_scaler = s_tree.array("P.pEL_CLEAN.scaler")
        P_pSTOF_scaler = s_tree.array("P.pSTOF.scaler")
        P_pPRHI_scaler = s_tree.array("P.PRHI.scaler")
        P_pPRLO_scaler = s_tree.array("P.PRLO.scaler")

        P_EDTM_scaler = s_tree.array("P.EDTM.scaler")
    
    NBCM = 5
    NTRIG = 6
    NPRE = 4
    NRATE = 6
    SHMSNRATE = 8

    bcm_name = ["BCM1 ", "BCM2 ", "BCM4A", "BCM4B", "BCM4C"]

    trig_name = ["TRIG1", "TRIG2", "TRIG3", "TRIG4", "TRIG5", "TRIG6"]

    PRE_name = ["40", "100", "150", "200"]

    rate_name = ["EL_LO_LO", "EL_LO", "EL_HI", "EL_REAL", "EL_CLEAN", "STOF"]

    SHMS_rate_name = ["EL_LO_LO", "EL_LO", "EL_HI",
                      "EL_REAL", "EL_CLEAN", "STOF", "PR_HI", "PR_LO"]

    bcm_value = [P_BCM1_scalerCharge, P_BCM2_scalerCharge,
                 P_BCM4A_scalerCharge, P_BCM4B_scalerCharge, P_BCM4C_scalerCharge]

    time_value = P_1Mhz_scalerTime

    current = [P_BCM1_scalerCurrent, P_BCM2_scalerCurrent,
                 P_BCM4A_scalerCurrent, P_BCM4B_scalerCurrent, P_BCM4C_scalerCurrent]

    trig_value = [P_pTRIG1_scaler, P_pTRIG2_scaler, P_pTRIG3_scaler,P_pTRIG4_scaler, P_pTRIG5_scaler, P_pTRIG6_scaler]
    #trig_value = [P_pTRIG1_scaler, P_pTRIG2_scaler, P_pEL_CLEAN_scaler,P_pTRIG4_scaler, P_pTRIG5_scaler, P_pTRIG6_scaler]

    acctrig_value = P_pL1ACCP_scaler

    PRE_value = [P_pPRE40_scaler, P_pPRE100_scaler,
                 P_pPRE150_scaler, P_pPRE200_scaler]

    SHMS_PRE_value = [P_pPRE40_scaler, P_pPRE100_scaler,
                      P_pPRE150_scaler, P_pPRE200_scaler]

    rate_value = [P_pEL_LO_LO_scaler, P_pEL_LO_scaler, P_pEL_HI_scaler,
                  P_pEL_REAL_scaler, P_pEL_CLEAN_scaler, P_pSTOF_scaler]

    SHMS_rate_value = [P_pEL_LO_LO_scaler, P_pEL_LO_scaler, P_pEL_HI_scaler,
                       P_pEL_REAL_scaler, P_pEL_CLEAN_scaler, P_pSTOF_scaler, P_pPRHI_scaler, P_pPRLO_scaler]

    EDTM_value = P_EDTM_scaler

    # Variables useful in Process
    # To find total charge
    name = [0]*NBCM
    charge_sum = [0]*NBCM
    time_sum = [0]*NBCM
    previous_charge = [0]*NBCM
    previous_time = 0
    previous_time = [0]*NBCM
    current_I = 0
    current_time = 0

    # To determine computer livetime
    name = [0]*NTRIG
    trig_sum = [0]*NTRIG
    previous_trig = [0]*NTRIG
    pretrigger = 0
    previous_pretrigger = 0
    acctrig_sum = 0
    previous_acctrig = 0

    # To determine HMS electronic livetime
    name = [0]*NPRE
    PRE_sum = [0]*NPRE
    previous_PRE = [0]*NPRE

    # To determine SHMS electronic livetime
    SHMS_PRE_sum = [0]*NPRE
    SHMS_previous_PRE = [0]*NPRE

    # To determine HMS trigger rates
    name = [0]*NRATE
    rate_sum = [0]*NRATE
    previous_rate = [0]*NRATE

    # To determine SHMS trigger rates
    rate_name = [0]*SHMSNRATE
    SHMS_rate_sum = [0]*SHMSNRATE
    SHMS_previous_rate = [0]*SHMSNRATE

    # To determine number of EDTM events
    EDTM_sum = 0
    EDTM_current = 0
    previous_EDTM = 0
    
    # Set bcm to use (0-bcm1, 1-bcm2, 2-bcm4A, 3-bcm4B, 4-bcm4C)
    bcm_ix = 1
    
    for ibcm in range(0, 5):
        previous_acctrig = (acctrig_value[0] - EDTM_current)
        previous_EDTM = EDTM_value[0]
        for itrig in range(0, NTRIG):
            previous_trig[itrig] = trig_value[itrig][0]
        for iPRE in range(0, NPRE):
            previous_PRE[iPRE] = PRE_value[iPRE][0]
            SHMS_previous_PRE[iPRE] = SHMS_PRE_value[iPRE][0]
        for iRATE in range(0, NRATE):
            previous_rate[iRATE] = rate_value[iRATE][0]
        for iRATE in range(0, SHMSNRATE):
            SHMS_previous_rate[iRATE] = SHMS_rate_value[iRATE][0]
        previous_time[ibcm] = time_value[0]
        previous_charge[ibcm] = bcm_value[ibcm][0]
        # Iterate over all scaler events to get various scaler values
        for i, evt in enumerate(s_evts):
            # Correction to bcm1 from Peter Bosted - removed this, b/c it's not for our data - NH 2024/04/08
            #if (current[ibcm][i] < 60):
            #    bcmcorr = 1.00+0.045*(math.log(60)-math.log(abs(current[ibcm][i]))/(math.log(60)-math.log(2)))
            #else:
            #    bcmcorr = 1.00+0.010*(current[ibcm][i]-60)/25
            #current[ibcm][i] = current[ibcm][i] * bcmcorr            
            if (time_value[i] != previous_time[ibcm]):
                # Current calculation using iterative charge and time values.
                # Iterate over current value then subtracting previous so that there is no double counting. Subtracted values are uncut.
                #current_I = (bcm_value[ibcm][i] - previous_charge[ibcm])/(time_value[i] - previous_time[ibcm])
                current_I = current[ibcm][i]
            if (abs( current[ibcm][i]-report_current) < thres_curr ):
                # Iterate over current value then subtracting previous so that there is no double counting. Subtracted values are uncut.
                charge_sum[ibcm] += (bcm_value[ibcm][i] - previous_charge[ibcm])
                time_sum[ibcm] += (time_value[i] - previous_time[ibcm])
            # Current cuts and selection of BCM1
            if (ibcm == bcm_ix and abs( current[ibcm][i]-report_current) < thres_curr):
                # EDTM scaler iteration.
                # Iterate over current value then subtracting previous so that there is no double counting. Subtracted values are uncut.
                EDTM_current = (EDTM_value[i] - previous_EDTM)
                EDTM_sum += EDTM_current
                # Accquired trigger sum calculation using iterative level 1 accepted values.
                # Iterate over current value then subtracting previous so that there is no double counting. Subtracted values are uncut.
                # Changing acctrig to not subtract EDTM to get CPULT for all events
                # acctrig_sum += ((acctrig_value[i] - EDTM_current) - previous_acctrig)
                acctrig_sum += ((acctrig_value[i]) - previous_acctrig)
                for itrig in range(0, NTRIG):
                    # Trigger scaler iteration.
                    # Iterate over current value then subtracting previous so that there is no double counting. Subtracted values are uncut.
                    trig_sum[itrig] += (trig_value[itrig][i] - previous_trig[itrig])
                    #print("trig_value[%s] = " %(itrig),trig_value[itrig][i])
                    #print("previous_trig[%s] = " %(itrig),previous_trig[itrig])
                for iPRE in range(0, NPRE):
                    # Pre-trig scaler iteration. Used in electronic LT calculations.
                    # Iterate over current value then subtracting previous so that there is no double counting. Subtracted values are uncut.
                    PRE_sum[iPRE] += (PRE_value[iPRE][i] - previous_PRE[iPRE])
                    SHMS_PRE_sum[iPRE] += (SHMS_PRE_value[iPRE][i] - SHMS_previous_PRE[iPRE])
                for iRATE in range(0, NRATE):
                    rate_sum[iRATE] += (rate_value[iRATE][i] - previous_rate[iRATE])
                for iRATE in range(0, SHMSNRATE):
                    SHMS_rate_sum[iRATE] += (SHMS_rate_value[iRATE][i] - SHMS_previous_rate[iRATE])
            # Changing acctrig to not subtract EDTM to get CPULT for all events
            # previous_acctrig = (acctrig_value[i] - EDTM_current)
            previous_acctrig = (acctrig_value[i])
            previous_EDTM = EDTM_value[i]
            for itrig in range(0, NTRIG):
                previous_trig[itrig] = trig_value[itrig][i]
            for iPRE in range(0, NPRE):
                previous_PRE[iPRE] = PRE_value[iPRE][i]
                SHMS_previous_PRE[iPRE] = SHMS_PRE_value[iPRE][i]
            for iRATE in range(0, NRATE):
                previous_rate[iRATE] = rate_value[iRATE][i]
            for iRATE in range(0, SHMSNRATE):
                SHMS_previous_rate[iRATE] = SHMS_rate_value[iRATE][i]
            previous_time[ibcm] = time_value[i]
            previous_charge[ibcm] = bcm_value[ibcm][i]

    # Define counter for trigger of interest
    for ps in PS_names:
        if ps is "PS1":
            shms_ps_ix = 0
        if ps is "PS2":
            shms_ps_ix = 1
        if ps is "PS3":
            hms_ps_ix = 2
        if ps is "PS4":
            hms_ps_ix = 3
        if ps is "PS5":
            coin_ps_ix = 4
        if ps is "PS6":
            coin_ps_ix = 5

    try:
        shms_ps_ix
    except NameError:
        shms_ps_ix = 0

    try:
        hms_ps_ix
    except NameError:
        hms_ps_ix = 2

    try:
        coin_ps_ix
    except NameError:
        coin_ps_ix = 4
        
    print("Debug: trig_sum", trig_sum)
    print("Debug: time: ", time_sum[bcm_ix], ", charge: ", charge_sum[bcm_ix], "\n\n")
    # Creates a dictionary for the calculated luminosity values 
    scalers = {
        "run number" : runNum,
        "time": time_sum[bcm_ix],
        "charge": charge_sum[bcm_ix],
        "curr_corr" : (charge_sum[bcm_ix]/time_sum[bcm_ix]-0.17)/(charge_sum[bcm_ix]/time_sum[bcm_ix]), # 0.17 uA current offset - NH 2024/05/04
        # "CPULT_scaler": acctrig_sum/((trig_sum[shms_ps_ix]/SHMS_PS) + (trig_sum[hms_ps_ix]/HMS_PS)), # GOOD
        #"CPULT_scaler": acctrig_sum/((trig_sum[shms_ps_ix]) + (trig_sum[hms_ps_ix]) - EDTM_sum),
        #"CPULT_scaler_uncern": (acctrig_sum/((trig_sum[shms_ps_ix]/SHMS_PS) + (trig_sum[hms_ps_ix]/HMS_PS)))*np.sqrt((1/(trig_sum[shms_ps_ix]/SHMS_PS))+(1/(trig_sum[hms_ps_ix]/HMS_PS))+(1/acctrig_sum)), # GOOD
        #"CPULT_scaler_uncern": (1/((trig_sum[shms_ps_ix]) + (trig_sum[hms_ps_ix])))*np.sqrt(acctrig_sum+EDTM_sum*2+((trig_sum[shms_ps_ix]) + (trig_sum[hms_ps_ix]))(acctrig_sum/((trig_sum[shms_ps_ix]) + (trig_sum[hms_ps_ix])))**2),
        #"CPULT_scaler_uncern": np.sqrt(((trig_sum[shms_ps_ix]) + (trig_sum[hms_ps_ix]))*.95*.05),
        "sent_edtm": EDTM_sum,
        "HMS_PS" : HMS_PS,
        "SHMS_PS" : SHMS_PS,
        "COIN_PS" : COIN_PS
            
    }

    if COIN_PS == None:
        if SHMS_PS == None:
            scalers.update({"sent_edtm_PS" : EDTM_sum/HMS_PS})
        elif HMS_PS == None:
            scalers.update({"sent_edtm_PS" : EDTM_sum/SHMS_PS})
        else:
            scalers.update({"sent_edtm_PS" : EDTM_sum/HMS_PS+EDTM_sum/SHMS_PS-EDTM_sum/(HMS_PS*SHMS_PS)})
    else:
        scalers.update({"sent_edtm_PS" : EDTM_sum/HMS_PS+EDTM_sum/SHMS_PS+EDTM_sum/COIN_PS+EDTM_sum/(SHMS_PS*HMS_PS*COIN_PS)-EDTM_sum/(HMS_PS*SHMS_PS)-EDTM_sum/(COIN_PS*SHMS_PS)-EDTM_sum/(HMS_PS*COIN_PS)})

    
    print("\n\n\n",trig_sum,"\n\n\n")

    if ("PS1" in PS_names or "PS2" in PS_names) and ("PS3" in PS_names or "PS4" in PS_names):
        print("Debug: in HMS + SHMS")
        scalers.update({"CPULT_scaler": acctrig_sum/((trig_sum[shms_ps_ix]/SHMS_PS) + (trig_sum[hms_ps_ix]/HMS_PS))})
        scalers.update({"CPULT_scaler_uncern": (acctrig_sum/((trig_sum[shms_ps_ix]/SHMS_PS) + (trig_sum[hms_ps_ix]/HMS_PS)))*np.sqrt((1/(trig_sum[shms_ps_ix]/SHMS_PS))+(1/(trig_sum[hms_ps_ix]/HMS_PS))+(1/acctrig_sum))})
    elif ("PS1" in PS_names or "PS2" in PS_names) and ("PS3" not in PS_names and "PS4" not in PS_names):
        print("Debug: in SHMS")
        scalers.update({"CPULT_scaler": acctrig_sum/((trig_sum[shms_ps_ix]/SHMS_PS))})
        scalers.update({"CPULT_scaler_uncern": (acctrig_sum/((trig_sum[shms_ps_ix]/SHMS_PS)))*np.sqrt((1/(trig_sum[shms_ps_ix]/SHMS_PS))+(1/acctrig_sum))})
    elif ("PS3" in PS_names or "PS4" in PS_names) and ("PS1" not in PS_names and "PS2" not in PS_names):
        print("Debug: in HMS")
        print("Debug: hms_ps_ix", hms_ps_ix, " trig_sum[hms_ps_ix] ", trig_sum[hms_ps_ix], " HMS_PS", HMS_PS)
        scalers.update({"CPULT_scaler": acctrig_sum/((trig_sum[2]/HMS_PS))})
        scalers.update({"CPULT_scaler_uncern": (acctrig_sum/((trig_sum[2]/HMS_PS)))*np.sqrt((1/(trig_sum[2]/HMS_PS))+(1/acctrig_sum))}) # this has sum bug that trig_sum[3] isn't being filled, don't know why - heinricn 2024/04/12
        scalers.update({})

    print("\nPre-scale values...")
    print("SHMS_PS : %s" % SHMS_PS, " HMS_PS : %s" % HMS_PS," COIN_PS : %s" % COIN_PS )
    scalers.update({"SHMSTRIG_scaler": trig_sum[shms_ps_ix]})
    scalers.update({"SHMS_eLT_scaler": 1 - ((6/5)*(SHMS_PRE_sum[1]-SHMS_PRE_sum[2])/(SHMS_PRE_sum[1]))})
    scalers.update({"SHMS_eLT_scaler_uncern": (SHMS_PRE_sum[1]-SHMS_PRE_sum[2])/(SHMS_PRE_sum[1])*np.sqrt((np.sqrt(SHMS_PRE_sum[1]) + np.sqrt(SHMS_PRE_sum[2]))/(SHMS_PRE_sum[1] - SHMS_PRE_sum[2]) + (np.sqrt(SHMS_PRE_sum[1])/SHMS_PRE_sum[1]))})
    scalers.update({"HMSTRIG_scaler": trig_sum[hms_ps_ix]})
    scalers.update({"HMS_eLT_scaler": 1 - ((6/5)*(PRE_sum[1]-PRE_sum[2])/(PRE_sum[1]))})
    scalers.update({"HMS_eLT_scaler_uncern": (PRE_sum[1]-PRE_sum[2])/(PRE_sum[1])*np.sqrt((np.sqrt(PRE_sum[1]) + np.sqrt(PRE_sum[2]))/(PRE_sum[1] - PRE_sum[2]) + (np.sqrt(PRE_sum[1])/PRE_sum[1]))})
    scalers.update({"COINTRIG_scaler": trig_sum[coin_ps_ix]})

    print("\n\nUsed current threshold value: %.2f uA" % thres_curr)

    for ibcm in range(0, NBCM):
        print("%s charge: %.3f uC, Beam over threshold for %.3f s" %
              (bcm_name[ibcm], charge_sum[ibcm], time_sum[ibcm]))

    print("\nL1ACC counts: %.0f, \nComputer Livetime: %f +/- %f" % (acctrig_sum, scalers["CPULT_scaler"],scalers["CPULT_scaler_uncern"]))
    for ps in PS_names:
        if ps == "PS1" or ps == "PS2":
            print("%s Prescaled Pretrigger Counts: %.0f" % (trig_name[0], scalers["SHMSTRIG_scaler"]))
            print("SHMS Electronic livetime: %f +/- %f" % (scalers["SHMS_eLT_scaler"], scalers["SHMS_eLT_scaler_uncern"]))            
        if ps == "PS3" or ps == "PS4":
            print("%s Prescaled Pretrigger Counts: %.0f" % (trig_name[0], scalers["HMSTRIG_scaler"]))
            print("HMS Electronic livetime: %f +/- %f" % (scalers["HMS_eLT_scaler"], scalers["HMS_eLT_scaler_uncern"]))
        if ps == "PS5" or ps == "PS6":
            print("%s Prescaled Pretrigger Counts: %.0f" % (trig_name[0], scalers["COINTRIG_scaler"]))

    print("EDTM Events: %.0f" % scalers["sent_edtm"])

    return scalers
