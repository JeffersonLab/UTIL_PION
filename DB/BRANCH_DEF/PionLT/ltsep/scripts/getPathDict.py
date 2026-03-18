#! /usr/bin/python

#
# Description: Script used in dynamic bash pathing. 
# I could simply just do this with aliases but this is 
# avoids any possible conflicts in naming conventions 
# for previously defined aliases
# ================================================================
# Time-stamp: "2022-06-30 02:25:22 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import sys
import ltsep as lt

CURRENT_ENV = sys.argv[1]

VOLATILEPATH=lt.SetPath(CURRENT_ENV).getPath("VOLATILEPATH")
ANALYSISPATH=lt.SetPath(CURRENT_ENV).getPath("ANALYSISPATH")
HCANAPATH=lt.SetPath(CURRENT_ENV).getPath("HCANAPATH")
REPLAYPATH=lt.SetPath(CURRENT_ENV).getPath("REPLAYPATH")
UTILPATH=lt.SetPath(CURRENT_ENV).getPath("UTILPATH")
PACKAGEPATH=lt.SetPath(CURRENT_ENV).getPath("PACKAGEPATH")
OUTPATH=lt.SetPath(CURRENT_ENV).getPath("OUTPATH")
ROOTPATH=lt.SetPath(CURRENT_ENV).getPath("ROOTPATH")
REPORTPATH=lt.SetPath(CURRENT_ENV).getPath("REPORTPATH")
CUTPATH=lt.SetPath(CURRENT_ENV).getPath("CUTPATH")
PARAMPATH=lt.SetPath(CURRENT_ENV).getPath("PARAMPATH")
SCRIPTPATH=lt.SetPath(CURRENT_ENV).getPath("SCRIPTPATH")
SIMCPATH=lt.SetPath(CURRENT_ENV).getPath("SIMCPATH")
LTANAPATH=lt.SetPath(CURRENT_ENV).getPath("LTANAPATH")
CACHEPATH=lt.SetPath(CURRENT_ENV).getPath("CACHEPATH")
ANATYPE=lt.SetPath(CURRENT_ENV).getPath("ANATYPE")
USER=lt.SetPath(CURRENT_ENV).getPath("USER")
HOST=lt.SetPath(CURRENT_ENV).getPath("HOST")

BashPathEntry=("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" % (VOLATILEPATH,ANALYSISPATH,HCANAPATH, REPLAYPATH, UTILPATH, PACKAGEPATH, OUTPATH, ROOTPATH, REPORTPATH, CUTPATH, PARAMPATH, SCRIPTPATH, ANATYPE, USER, HOST, SIMCPATH, LTANAPATH, CACHEPATH))
print(BashPathEntry)
