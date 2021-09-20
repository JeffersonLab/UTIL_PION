#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2021-08-31 00:27:50 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import uproot as up
import numpy as np


def scaler(PS_names, SHMS_PS, HMS_PS, thres_curr,report_current,REPLAYPATH,runNum,MaxEvent,s_tree,s_branch):

    '''
    SCALER TREE, TSP
    '''

    # s_evts = len(s_tree)
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

    trig_value = [P_pTRIG1_scaler, P_pTRIG2_scaler, P_pTRIG3_scaler,
                  P_pTRIG4_scaler, P_pTRIG5_scaler, P_pTRIG6_scaler]

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
        for i, evt in enumerate(s_evts):
            if (time_value[i] != previous_time[ibcm]):
                current_I = (bcm_value[ibcm][i] -
                             previous_charge[ibcm])/(time_value[i] - previous_time[ibcm])
            if (abs( current[ibcm][i]-report_current) < thres_curr ):
                charge_sum[ibcm] += (bcm_value[ibcm][i] - previous_charge[ibcm])
                time_sum[ibcm] += (time_value[i] - previous_time[ibcm])
            if (ibcm == 2 and abs( current[ibcm][i]-report_current) < thres_curr):
                EDTM_current = (EDTM_value[i] - previous_EDTM)
                EDTM_sum += EDTM_current
                acctrig_sum += ((acctrig_value[i] - EDTM_current) - previous_acctrig)
                for itrig in range(0, NTRIG):
                    trig_sum[itrig] += (trig_value[itrig][i] - previous_trig[itrig])
                    # print("trig_value[%s] = " %(itrig),trig_value[itrig][i])
                for iPRE in range(0, NPRE):
                    PRE_sum[iPRE] += (PRE_value[iPRE][i] - previous_PRE[iPRE])
                    SHMS_PRE_sum[iPRE] += (SHMS_PRE_value[iPRE][i] - SHMS_previous_PRE[iPRE])
                for iRATE in range(0, NRATE):
                    rate_sum[iRATE] += (rate_value[iRATE][i] - previous_rate[iRATE])
                for iRATE in range(0, SHMSNRATE):
                    SHMS_rate_sum[iRATE] += (SHMS_rate_value[iRATE][i] - SHMS_previous_rate[iRATE])
            previous_acctrig = (acctrig_value[i] - EDTM_current)
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

    if PS_names[0] is "PS1":
        shms_ps_ix = 0
    if PS_names[0] is "PS2":
        shms_ps_ix = 1
    if PS_names[1] is "PS3":
        hms_ps_ix = 2
    if PS_names[1] is "PS4":
        hms_ps_ix = 3
        
    scalers = {
        "run number" : runNum,
        "%s" % PS_names[shms_ps_ix]: SHMS_PS,
        "%s" % PS_names[1]: HMS_PS,
        "time": time_sum[2],
        "charge": charge_sum[2],
        "SHMSTRIG_scaler": trig_sum[shms_ps_ix],
        "HMSTRIG_scaler": trig_sum[hms_ps_ix],
#        "CPULT_scaler": acctrig_sum/((trig_sum[shms_ps_ix]/SHMS_PS) + (trig_sum[hms_ps_ix]/HMS_PS)),
        "CPULT_scaler": 1-acctrig_sum/((trig_sum[shms_ps_ix]) + (trig_sum[hms_ps_ix])),
        "CPULT_scaler_uncern": (acctrig_sum/((trig_sum[shms_ps_ix]/SHMS_PS) + (trig_sum[hms_ps_ix]/HMS_PS)))*np.sqrt((1/(trig_sum[shms_ps_ix]/SHMS_PS))+(1/(trig_sum[hms_ps_ix]/HMS_PS))+(1/acctrig_sum)),
        "HMS_eLT": 1 - ((6/5)*(PRE_sum[1]-PRE_sum[2])/(PRE_sum[1])),
        "HMS_eLT_uncern": (PRE_sum[1]-PRE_sum[2])/(PRE_sum[1])*np.sqrt((np.sqrt(PRE_sum[1]) + np.sqrt(PRE_sum[2]))/(PRE_sum[1] - PRE_sum[2]) + (np.sqrt(PRE_sum[1])/PRE_sum[1])),
        "SHMS_eLT": 1 - ((6/5)*(SHMS_PRE_sum[1]-SHMS_PRE_sum[2])/(SHMS_PRE_sum[1])),
        "SHMS_eLT_uncern": (SHMS_PRE_sum[1]-SHMS_PRE_sum[2])/(SHMS_PRE_sum[1])*np.sqrt((np.sqrt(SHMS_PRE_sum[1]) + np.sqrt(SHMS_PRE_sum[2]))/(SHMS_PRE_sum[1] - SHMS_PRE_sum[2]) + (np.sqrt(SHMS_PRE_sum[1])/SHMS_PRE_sum[1])),
        "sent_edtm": EDTM_sum
            
    }


    print("Using prescale factors: %s %.0f, %s %.0f\n" % (PS_names[0],SHMS_PS,PS_names[1],HMS_PS))
    print("\n\nUsed current threshold value: %.2f uA" % thres_curr)

    for ibcm in range(0, NBCM):
        print("%s charge: %.3f uC, Beam over threshold for %.3f s" %
              (bcm_name[ibcm], charge_sum[ibcm], time_sum[ibcm]))

    print("\n\n")

    print("L1ACC counts: %.0f, \n%s Prescaled Pretrigger Counts: %.0f \n%s Prescaled Pretrigger Counts: %.0f \nComputer Livetime: %f +/- %f" %
          (acctrig_sum, trig_name[0], scalers["SHMSTRIG_scaler"], trig_name[2], scalers["HMSTRIG_scaler"], scalers["CPULT_scaler"], scalers["CPULT_scaler_uncern"]))

    print("HMS Electronic livetime: %f +/- %f" %
          (scalers["HMS_eLT"], scalers["HMS_eLT_uncern"]))

    print("SHMS Electronic livetime: %f +/- %f" %
          (scalers["SHMS_eLT"], scalers["SHMS_eLT_uncern"]))

    print("EDTM Events: %.0f" % scalers["sent_edtm"])

    return scalers
