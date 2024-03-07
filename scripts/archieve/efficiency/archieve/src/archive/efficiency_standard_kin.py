#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2022-09-08 06:13:02 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

################################################################################################################################################
'''
'''

def dictionary(ANATYPE,REPLAYPATH,runNum):

    # Open report file to grab prescale values and tracking efficiency
    standard_kinematics = REPLAYPATH+"/DBASE/COIN/DB_KaonLT/standard_Offline%sLT.kinematics" % ANATYPE

    with open(standard_kinematics) as f:
        effDict = {

            'gpbeam': None,
            'hpcentral': None,
            'htheta_lab': None,
            'hpartmass': None,
            'ppcentral': None,
            'ptheta_lab': None,
            'ppartmass': None,

        }

        runcheck = False
        # Search for keywords, then save as value in dictionary
        for line in f:
            if "#" in line:
                continue
            else:
                if '-' in line and '=' not in line:
                    runline = line.split('-')
                    #print(runline)
                    if int(runNum) >= int(runline[0]) and int(runNum) <= int(runline[1]):
                        runcheck = True
                    else:
                        runcheck = False
                if runcheck:
                        data = line.split('=')
                        for key,val in effDict.items():
                            #print(data[0])
                            if key in data[0]:
                                effDict[key] = float(data[1])
                    
        print(effDict)

    return effDict
