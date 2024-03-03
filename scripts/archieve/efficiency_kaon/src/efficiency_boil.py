#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-12-18 16:42:12 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import math

def dictionary(BCM_current,DEBUG=False):

    I = BCM_current
    dI = 0.1e-3 # 0.1 uA (https://hallcweb.jlab.org/DocDB/0010/001034/002/BCM1and2StabilityFall18toSummer19v3.pdf)
    m0 = -7.899e-4
    dm0 = 1.829e-4
    
    boil_eff = 1 - abs(m0)*I
    boil_eff_err = math.sqrt((I**2)*(dm0**2)+(m0**2)*(dI**2))
    
    effDict = {

        "BOIL_Eff" : boil_eff,
        "BOIL_Eff_ERROR" : boil_eff_err
    }
    
    return effDict
