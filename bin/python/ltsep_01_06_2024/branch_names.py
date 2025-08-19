#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-01-08 18:30:07 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

# Dictionary of all used branches
# Add more if required (note: make sure to add to DB/BRANCH_DEF/<ANATYPE>LT/<RUNTYPE_FILE> also)
branch_dict = {
    # HMS info
    "H_dc_InsideDipoleExit" : "H.dc.InsideDipoleExit",
    "H_hod_goodscinhit" : "H.hod.goodscinhit",
    "H_hod_goodstarttime" : "H.hod.goodstarttime",
    # Beta is velocity of particle between pairs of hodoscopes
    "H_gtr_beta" : "H.gtr.beta",
    "H_dc_x_fp" : "H.dc.x_fp",
    "H_dc_y_fp" : "H.dc.y_fp",
    "H_dc_xp_fp" : "H.dc.xp_fp",
    "H_dc_yp_fp" : "H.dc.yp_fp",
    # xpfp -> Theta
    "H_gtr_xp" : "H.gtr.th",
    # ypfp -> Phi
    "H_gtr_yp" : "H.gtr.ph",
    # dp is Delta
    "H_gtr_dp" : "H.gtr.dp",
    "H_gtr_p" : "H.gtr.p",
    "H_cal_etotnorm" : "H.cal.etotnorm",
    "H_cal_etottracknorm" : "H.cal.etottracknorm",
    "H_cer_npeSum" : "H.cer.npeSum",
    "H_W" : "H.kin.primary.W",
    "H_cal_etotnorm" : "H.cal.etotnorm",
    "H_cer_npeSum" : "H.cer.npeSum",
    "H_gtr_dp" : "H.gtr.dp",
    "H_tr_tg_th" : "H.gtr.th",
    "H_tr_tg_ph" : "H.gtr.ph",
    "H_gtr_beta" : "H.gtr.beta",
    "H_tr_chi2" : "H.tr.chi2",
    "H_tr_ndof" : "H.tr.ndof",
    "H_hod_goodscinhit" : "H.hod.goodscinhit",
    "H_hod_betanotrack" : "H.hod.betanotrack",
    "H_hod_goodstarttime" : "H.hod.goodstarttime",
    "H_dc_ntrack" : "H.dc.ntrack",
    "H_dc_1x1_nhit" : "H.dc.1x1.nhit",
    "H_dc_1u2_nhit" : "H.dc.1u2.nhit",
    "H_dc_1u1_nhit" : "H.dc.1u1.nhit",
    "H_dc_1v1_nhit" : "H.dc.1v1.nhit",
    "H_dc_1x2_nhit" : "H.dc.1x2.nhit",
    "H_dc_1v2_nhit" : "H.dc.1v2.nhit",
    "H_dc_2x1_nhit" : "H.dc.2x1.nhit",
    "H_dc_2u2_nhit" : "H.dc.2u2.nhit",
    "H_dc_2u1_nhit" : "H.dc.2u1.nhit",
    "H_dc_2v1_nhit" : "H.dc.2v1.nhit",
    "H_dc_2x2_nhit" : "H.dc.2x2.nhit",
    "H_dc_2v2_nhit" : "H.dc.2v2.nhit",

    # SHMS info
    "P_cal_fly_earray" : "P.cal.fly.earray",
    "P_cal_pr_eplane" : "P.cal.pr.eplane",
    "P_cal_etotnorm" : "P.cal.etotnorm",
    "P_aero_npeSum" : "P.aero.npeSum",
    "P_hgcer_npeSum" : "P.hgcer.npeSum",
    "P_hgcer_xAtCer" : "P.hgcer.xAtCer",
    "P_hgcer_yAtCer" : "P.hgcer.yAtCer",
    "P_ngcer_npeSum" : "P.ngcer.npeSum",
    "P_ngcer_xAtCer" : "P.ngcer.xAtCer",
    "P_ngcer_yAtCer" : "P.ngcer.yAtCer",
    "P_aero_xAtCer" : "P.aero.xAtAero",
    "P_aero_yAtCer" : "P.aero.yAtAero",
    "P_dc_InsideDipoleExit" : "P.dc.InsideDipoleExit",
    "P_hod_goodscinhit" :     "P.hod.goodscinhit",
    "P_hod_goodstarttime" : "P.hod.goodstarttime",
    # Beta is velocity of particle between pairs of hodoscopes
    "P_gtr_beta" : "P.gtr.beta",
    "P_gtr_x" : "P.gtr.x",
    "P_gtr_y" : "P.gtr.y",
    "P_dc_x_fp" : "P.dc.x_fp",
    "P_dc_y_fp" : "P.dc.y_fp",
    "P_dc_xp_fp" : "P.dc.xp_fp",
    "P_dc_yp_fp" : "P.dc.yp_fp",
    # xpfp -> Theta
    "P_gtr_xp" : "P.gtr.th",
    # ypfp -> Phi
    "P_gtr_yp" : "P.gtr.ph",
    "P_gtr_p" : "P.gtr.p",
    # dp is Delta
    "P_gtr_dp" : "P.gtr.dp",
    "P_cal_etotnorm" : "P.cal.etotnorm",
    "P_cal_etottracknorm" : "P.cal.etottracknorm",
    "P_aero_npeSum" : "P.aero.npeSum",
    "P_aero_xAtAero" : "P.aero.xAtAero",
    "P_aero_yAtAero" : "P.aero.yAtAero",
    "P_hgcer_npeSum" : "P.hgcer.npeSum",
    "P_hgcer_xAtCer" : "P.hgcer.xAtCer",
    "P_hgcer_yAtCer" : "P.hgcer.yAtCer",
    "P_ngcer_npeSum" : "P.ngcer.npeSum",
    "P_ngcer_xAtCer" : "P.ngcer.xAtCer",
    "P_ngcer_yAtCer" : "P.ngcer.yAtCer",
    "P_cal_etotnorm" : "P.cal.etotnorm",
    "P_hgcer_npeSum" : "P.hgcer.npeSum",
    "P_ngcer_npeSum" : "P.ngcer.npeSum",
    "P_aero_npeSum" : "P.aero.npeSum",
    "P_gtr_dp" : "P.gtr.dp",
    "P_gtr_th" : "P.gtr.th",
    "P_gtr_ph" : "P.gtr.ph",
    "P_gtr_beta" : "P.gtr.beta",
    "P_tr_chi2" : "P.tr.chi2",
    "P_tr_ndof" : "P.tr.ndof",
    "P_hod_goodscinhit" : "P.hod.goodscinhit",
    "P_hod_betanotrack" : "P.hod.betanotrack",
    "P_hod_goodstarttime" : "P.hod.goodstarttime",
    "P_dc_ntrack" : "P.dc.ntrack",
    "P_dc_1x1_nhit" : "P.dc.1x1.nhit",
    "P_dc_1u2_nhit" : "P.dc.1u2.nhit",
    "P_dc_1u1_nhit" : "P.dc.1u1.nhit",
    "P_dc_1v1_nhit" : "P.dc.1v1.nhit",
    "P_dc_1x2_nhit" : "P.dc.1x2.nhit",
    "P_dc_1v2_nhit" : "P.dc.1v2.nhit",
    "P_dc_2x1_nhit" : "P.dc.2x1.nhit",
    "P_dc_2u2_nhit" : "P.dc.2u2.nhit",
    "P_dc_2u1_nhit" : "P.dc.2u1.nhit",
    "P_dc_2v1_nhit" : "P.dc.2v1.nhit",
    "P_dc_2x2_nhit" : "P.dc.2x2.nhit",
    "P_dc_2v2_nhit" : "P.dc.2v2.nhit",

    # Raster
    "raster_x" : "P.rb.x",
    "raster_y" : "P.rb.y",
    "raster_z" : "P.rb.z",

    # BPM target
    "bpm_tar_x" : "P.rb.raster.fr_xbpm_tar",
    "bpm_tar_y" : "P.rb.raster.fr_ybpm_tar",

    # Kinematic quantitites
    "Q2" : "H.kin.primary.Q2",
    "W" : "H.kin.primary.W",
    "epsilon" : "H.kin.primary.epsilon",
    "ph_q" : "P.kin.secondary.ph_xq",
    "ph_recoil" : "P.kin.secondary.ph_bq",
    "th_q" : "P.kin.secondary.th_xq",
    "th_recoil" : "P.kin.secondary.th_bq",
    "emiss" : "P.kin.secondary.emiss",
    "MMpi" : "P.kin.secondary.MMpi",
    "MMK" : "P.kin.secondary.MMK",
    "MMp" : "P.kin.secondary.MMp",
    "MandelT" : "P.kin.secondary.MandelT",
    "MandelU" : "P.kin.secondary.MandelU",
    "pmiss" : "P.kin.secondary.pmiss",
    "pmiss_x" : "P.kin.secondary.pmiss_x",
    "pmiss_y" : "P.kin.secondary.pmiss_y",
    "pmiss_z" : "P.kin.secondary.pmiss_z",
    "Erecoil" : "P.kin.secondary.Erecoil",
    "emiss_nuc" : "P.kin.secondary.emiss_nuc",
    "Mrecoil" : "P.kin.secondary.Mrecoil",

    # Current
    "H_bcm_bcm1_AvgCurrent" : "H.bcm.bcm1.AvgCurrent",
    "H_bcm_bcm2_AvgCurrent" : "H.bcm.bcm2.AvgCurrent",
    "H_bcm_bcm4a_AvgCurrent" : "H.bcm.bcm4a.AvgCurrent",
    "H_bcm_bcm4b_AvgCurrent" : "H.bcm.bcm4b.AvgCurrent",
    "H_bcm_bcm4c_AvgCurrent" : "H.bcm.bcm4c.AvgCurrent",

    # Timing info
    "CTime_eKCoinTime_ROC1" : "CTime.eKCoinTime_ROC1",
    "CTime_ePiCoinTime_ROC1" : "CTime.ePiCoinTime_ROC1",
    "CTime_epCoinTime_ROC1" : "CTime.epCoinTime_ROC1",

    "P_RF_tdcTime" : "T.coin.pRF_tdcTime",
    "P_hod_fpHitsTime" : "P.hod.fpHitsTime",
    "H_RF_Dist" : "RFTime.HMS_RFtimeDist",
    "P_RF_Dist" : "RFTime.SHMS_RFtimeDist",


    "T_coin_pTRIG1_ROC1_tdcTimeRaw" : "T.coin.pTRIG1_ROC1_tdcTimeRaw",
    "T_coin_pTRIG1_ROC2_tdcTimeRaw" : "T.coin.pTRIG1_ROC2_tdcTimeRaw",
    "T_coin_pTRIG1_ROC1_tdcTime" : "T.coin.pTRIG1_ROC1_tdcTime",
    "T_coin_pTRIG1_ROC2_tdcTime" : "T.coin.pTRIG1_ROC2_tdcTime",

    "T_coin_pTRIG2_ROC1_tdcTimeRaw" : "T.coin.pTRIG2_ROC1_tdcTimeRaw",
    "T_coin_pTRIG2_ROC2_tdcTimeRaw" : "T.coin.pTRIG2_ROC2_tdcTimeRaw",
    "T_coin_pTRIG2_ROC1_tdcTime" : "T.coin.pTRIG2_ROC1_tdcTime",
    "T_coin_pTRIG2_ROC2_tdcTime" : "T.coin.pTRIG2_ROC2_tdcTime",

    "T_coin_pTRIG3_ROC1_tdcTimeRaw" : "T.coin.pTRIG3_ROC1_tdcTimeRaw",
    "T_coin_pTRIG3_ROC2_tdcTimeRaw" : "T.coin.pTRIG3_ROC2_tdcTimeRaw",
    "T_coin_pTRIG3_ROC1_tdcTime" : "T.coin.pTRIG3_ROC1_tdcTime",
    "T_coin_pTRIG3_ROC2_tdcTime" : "T.coin.pTRIG3_ROC2_tdcTime",

    "T_coin_pTRIG4_ROC1_tdcTimeRaw" : "T.coin.pTRIG4_ROC1_tdcTimeRaw",
    "T_coin_pTRIG4_ROC2_tdcTimeRaw" : "T.coin.pTRIG4_ROC2_tdcTimeRaw",
    "T_coin_pTRIG4_ROC1_tdcTime" : "T.coin.pTRIG4_ROC1_tdcTime",
    "T_coin_pTRIG4_ROC2_tdcTime" : "T.coin.pTRIG4_ROC2_tdcTime",

    "T_coin_pTRIG5_ROC1_tdcTimeRaw" : "T.coin.pTRIG5_ROC1_tdcTimeRaw",
    "T_coin_pTRIG5_ROC2_tdcTimeRaw" : "T.coin.pTRIG5_ROC2_tdcTimeRaw",
    "T_coin_pTRIG5_ROC1_tdcTime" : "T.coin.pTRIG5_ROC1_tdcTime",
    "T_coin_pTRIG5_ROC2_tdcTime" : "T.coin.pTRIG5_ROC2_tdcTime",

    "T_coin_pTRIG6_ROC1_tdcTimeRaw" : "T.coin.pTRIG6_ROC1_tdcTimeRaw",
    "T_coin_pTRIG6_ROC2_tdcTimeRaw" : "T.coin.pTRIG6_ROC2_tdcTimeRaw",
    "T_coin_pTRIG6_ROC1_tdcTime" : "T.coin.pTRIG6_ROC1_tdcTime",
    "T_coin_pTRIG6_ROC2_tdcTime" : "T.coin.pTRIG6_ROC2_tdcTime",

    "T_coin_pFADC_TREF_ROC2_adcPed" : "T.coin.pFADC_TREF_ROC2_adcPed",
    "T_coin_hFADC_TREF_ROC1_adcPed" : "T.coin.hFADC_TREF_ROC1_adcPed",
    "T_coin_pFADC_TREF_ROC2_adcPulseTimeRaw" : "T.coin.pFADC_TREF_ROC2_adcPulseTimeRaw",
    "T_coin_hFADC_TREF_ROC1_adcPulseTimeRaw" : "T.coin.hFADC_TREF_ROC1_adcPulseTimeRaw",
    "T_coin_pEDTM_tdcTimeRaw" : "T.coin.pEDTM_tdcTimeRaw",
    "T_coin_pEDTM_tdcTime" : "T.coin.pEDTM_tdcTime",

    # Misc quantities
    "RFFreq" : "MOFC1FREQ",
    "RFFreqDiff" : "MOFC1DELTA",
    "EvtType" : "fEvtHdr.fEvtType",
}
