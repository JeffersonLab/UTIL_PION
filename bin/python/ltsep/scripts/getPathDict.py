#! /usr/bin/python

#
# Description: 
# ================================================================
# Time-stamp: "2021-11-03 05:52:22 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import sys
import ltsep as lt

CURRENT_ENV = sys.argv[1]

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
ANATYPE=lt.SetPath(CURRENT_ENV).getPath("ANATYPE")
USER=lt.SetPath(CURRENT_ENV).getPath("USER")
HOST=lt.SetPath(CURRENT_ENV).getPath("HOST")

BashPathEntry=("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" % (HCANAPATH, REPLAYPATH, UTILPATH, PACKAGEPATH, OUTPATH, ROOTPATH, REPORTPATH, CUTPATH, PARAMPATH, SCRIPTPATH, ANATYPE, USER, HOST))
print(BashPathEntry)
