#! /usr/bin/python

#
# Description: This package will perform many tasks required for physics analysis in hall c
#              (although can also be expanded to any analysis using root files to varying success)
# Analysis script required format for applying cuts...

'''
import uproot as up
sys.path.insert(0, 'path_to/bin/python/')
import kaonlt as klt

# Convert root leaf to array with uproot
# Array name must match what is defined in DB/CUTS/general/
array  = tree.array("leaf")

# Not required for applying cuts, but required for converting back to root files
r = klt.pyRoot()

fout = "<path_to_run_type_cut>"
c = klt.pyPlot(None) # See below for pyPlot class definition
readDict = c.read_dict(fout) # read in run type cuts file and makes dictionary

# This method calls several methods in kaonlt package. It is required to create properly formated
# dictionaries. The evaluation must be in the analysis script because the analysis variables (i.e. the
# leaves of interest) are not defined in the kaonlt package. This makes the system more flexible
# overall, but a bit more cumbersome in the analysis script. Perhaps one day a better solution will be
# implimented.
def make_cutDict(cut,inputDict=None):

    global c

    c = klt.pyPlot(readDict)
    x = c.w_dict(cut)
    
    # Only for first key of dictionary
    if inputDict == None:
        inputDict = {}
        
    # Update dictionary with cuts (as strings) from readDict
    for key,val in readDict.items():
        if key == cut:
            inputDict.update({key : {}})

    # Evaluate strings to cut values. Creates a dictionary in a dictionary...dict-ception!
    for i,val in enumerate(x):
        tmp = x[i]
        # Checks for removed leaves
        if tmp == "":
            continue
        else:
            inputDict[cut].update(eval(tmp))
        
    return inputDict

cutDict = make_cutDict("cut1")
cutDict = make_cutDict("cut2",cutDict)
# Continue this for all run type cuts required

# ---> If multple run type files are required then define a new run type file altogether. Do not try to 
# chain run type files. It can be done, but is computationally wasteful and pointless.

# To apply cuts to array...
c.add_cut(array,"cut#")

'''
# ================================================================
# Time-stamp: "2020-05-02 15:00:37 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#


from __future__ import division
import logging

# Gets rid of matplot logging DEBUG messages
plt_logger = logging.getLogger('matplotlib')
plt_logger.setLevel(logging.WARNING)

# Suppresses unwanted numpy warning
import warnings
import numpy as np
warnings.simplefilter(action='ignore', category=FutureWarning)

from ROOT import TFile, TH1F
import matplotlib.pyplot as plt
from matplotlib import interactive
from matplotlib import colors
import uproot as up
import pandas as pd
from csv import DictReader
import time, math, sys, subprocess

# garbage collector
import gc
gc.collect()

'''
When calling kaonlt package, you may define a dictionary in the script. This dictionary will contain
the cuts of interest (defined in a CUTS directory).  These cuts are read in through the read_dict()
method of the pyPlot() class. The pyDict class is not explicitly called, but rather called implicitly
by other classes.
'''
class pyDict(dict):
    
    def __init__(self,inputTree):
        self.inputTree = inputTree
        
'''
This class, with its findBranch method, will grab the leaves in a branch using uproot package. This
takes the tree as an input.
'''
class pyBranch(pyDict):

    def findBranch(self,inputBranch, inputLeaf):
        tree = self.inputTree
        branch = tree.array(inputBranch)
        branch  = list(zip(*branch)) # Match elements to proper leaves
        leafList = tree[inputBranch].interpretation.fromdtype.descr
        i=0
        for name,typ in leafList:
            if name == inputLeaf:
                leaf = name
                leafVal = i
                break
            i+=1
        leafHist = branch[leafVal]

        return np.array(leafHist)

'''    
This class is for converting files into root files after the analysis steps
'''
class pyRoot():

    # Save arrays,lists,etc. from csv to root file as histograms
    def csv2root(self,inputDict,rootName):
        try:
            tmp = ""
            hist_key = []*len(inputDict)
            hist_val = []*len(inputDict)
            for key,val in inputDict.items():
                tmp = "hist_%s" % key
                tmp = TH1F( tmp, '%s' % key, len(val), 0., max(val))
                hist_key.append(tmp)
                hist_val.append(val)

            f = TFile( rootName, 'recreate' )
            for i, evt in enumerate(hist_val):
                for j, hevt in enumerate(hist_val[i]):
                    print(hist_key[i], "-> ", hevt)
                    hist_key[i].Fill(hevt)
                hist_key[i].Write()
 
            f.Write()
            f.Close()
        except TypeError:
            print("\nERROR 1: Only current accepting 1D array/list values\n")

'''            
This class stores a variety of equations often used in the KaonLT analysis procedure
'''
class pyEquation():

    # Define missing mass calculation
    def missmass():
        print("missmass")

'''
This is the most extensive class of the kaonlt package. This class will perform many required tasks
for doing in depth analysis in python. This class does not require, but will use the pyDict class to
apply cuts. Set the dictionary to None if no cuts are required.
'''
class pyPlot(pyDict):
    
    def __init__(self, REPLAYPATH,cutDict=None,DEBUG=False):
        self.REPLAYPATH = REPLAYPATH
        self.cutDict = cutDict
        self.DEBUG = DEBUG

    # A method for defining a bin. This may be called in any matplotlib package plots.
    # This will calculate a suitable bin width and use that to equally distribute the bin size
    def setbin(self,plot,numbin,xmin=None,xmax=None):
        
        if (xmin or xmax):
            leaf = self.fixBin(plot,plot,xmin,xmax)
        else:
            leaf = plot
            
        binwidth = (abs(leaf).max()-abs(leaf).min())/numbin
        
        bins = np.arange(min(leaf), max(leaf) + binwidth, binwidth)

        return bins

    # This method is complimentary to setbin(). This will cut the distribution based off the min and max
    # array values
    def fixBin(self,cut,plot,low,high):
            
        arrCut = cut
        arrPlot = plot
        arrPlot = arrPlot[(arrCut > low) & (arrCut < high)]

        return arrPlot

    def cut_RF(self,runNum,MaxEvent):
        TimingCutFile = self.REPLAYPATH+'/UTIL_PION/DB/PARAM/Timing_Parameters.csv'
        # rootName = "/lustre19/expphy/volatile/hallc/c-pionlt/sjdkay/ROOTfiles/Proton_Analysis/Pass3/Proton_coin_replay_production_%s_%s.root" % (self.REPLAYPATH, runNum, MaxEvent)
        rootName = "%s/UTIL_PION/ROOTfiles/coin_replay_Full_Lumi_%s_%s.root" % (self.REPLAYPATH,runNum,MaxEvent)
        e_tree = up.open(rootName)["T"]
        TimingCutf = open(TimingCutFile)
        PromptPeak = [0, 0, 0]
        linenum = 0 # Count line number we're on
        TempPar = -1 # To check later
        for line in TimingCutf: # Read all lines in the cut file
            linenum += 1 # Add one to line number at start of loop
            if(linenum > 1): # Skip first line
                line = line.partition('#')[0] # Treat anything after a # as a comment and ignore it
                line = line.rstrip()
                array = line.split(",") # Convert line into an array, anything after a comma is a new entry
                if(int(runNum) in range (int(array[0]), int(array[1])+1)): # Check if run number for file is within any of the ranges specified in the cut file
                    TempPar += 2 # If run number is in range, set to non -1 value
                    BunchSpacing = float(array[2]) # Bunch spacing in ns
                    RF_Offset = float(array[9]) # Offset for RF timing cut
        TimingCutf.close() # After scanning all lines in file, close file
        if(TempPar == -1): # If value is still -1, run number provided didn't match any ranges specified so exit
            print("!!!!! ERROR !!!!!\n Run number specified does not fall within a set of runs for which cuts are defined in %s\n!!!!! ERROR !!!!!" % TimingCutFile)
            sys.exit(3)
        elif(TempPar > 1):
            print("!!! WARNING!!! Run number was found within the range of two (or more) line entries of %s !!! WARNING !!!" % TimingCutFile)
            print("The last matching entry will be treated as the input, you should ensure this is what you want")
        P_RF_tdcTime = e_tree.array("T.coin.pRF_tdcTime")
        P_hod_fpHitsTime = e_tree.array("P.hod.fpHitsTime")
        RF_CutDist = np.array([ ((RFTime-StartTime + RF_Offset)%(BunchSpacing)) for (RFTime, StartTime) in zip(P_RF_tdcTime, P_hod_fpHitsTime)]) # In python x % y is taking the modulo y of x

    # This method reads in the CUTS and converts them to a dictionary. 
    def read_dict(self,fout,runNum):

        # Open run type cuts of interest
        f = open(fout)
        cutDict = {}
        for line in f:
            if "#" in line:
                continue
            else:
                
                line = line.split("=",1)
                # Grab run type cut name
                typName = line[0].rstrip()
                typName = typName.lstrip()
                # Grab run type cuts required, note at this stage the cuts to be removed are bunched
                # together still
                typCuts = line[1].split("+")
                if (self.DEBUG):
                    print("Type ", typName)
                    print("Cuts ", typCuts)
                
                # Loop over run type cuts
                for i,evt in enumerate(typCuts):
                    # Split any cuts to be removed
                    minusCuts = evt.split("-")
                    # Ignore first element, since it will always be an added cut
                    minus = minusCuts[1:]
                    # Define first cut to be added, any other cuts to be added will be done in future
                    # iteration over run type cuts
                    cutplus = minusCuts[0].rstrip()
                    cutplus = cutplus.lstrip()
                    if (self.DEBUG):
                        print("+ ",cutplus)

                    ##############
                    # Added cuts #
                    ##############
                    
                    # Matches run type cuts with the general cuts (e.g pid, track, etc.)
                    if "pid" in cutplus:
                        plusfout = self.REPLAYPATH+"/UTIL_PION/DB/CUTS/general/pid.cuts"
                    elif "track" in cutplus:
                        plusfout = self.REPLAYPATH+"/UTIL_PION/DB/CUTS/general/track.cuts"
                    elif "accept" in cutplus:
                        plusfout = self.REPLAYPATH+"/UTIL_PION/DB/CUTS/general/accept.cuts"
                    elif "coin_time" in cutplus:
                        plusfout = self.REPLAYPATH+"/UTIL_PION/DB/CUTS/general/coin_time.cuts"
                    elif "current" in cutplus:
                        plusfout = self.REPLAYPATH+"/UTIL_PION/DB/CUTS/general/current.cuts"
                    elif "misc" in cutplus:
                        plusfout = self.REPLAYPATH+"/UTIL_PION/DB/CUTS/general/misc.cuts"
                    else:
                        print("!!!!ERROR!!!!: Added cut %s not defined in /UTIL_PION/DB/CUTS/general/" % cutplus) # ERROR 2
                        print("Cut must be pid, track, accept, coin_time or current")
                        continue
                    cutplus = cutplus.split(".")
                    if len(cutplus) == 2:
                        cutplus = str(cutplus[1])
                        # print("cutplus ", cutplus)
                    elif len(cutplus) > 2:
                        cutplus = str(cutplus[2])
                        # print("cutplus ", cutplus)
                    else:
                        # print("ERROR 5: %s cut not found in %s" % (cutplus,plusfout))
                        continue

                    # Open general cuts file of interest to be added to dictionary
                    fplus = open(plusfout)
                    for lplus in fplus:
                        if "#" in lplus:
                            continue
                        else:
                            lplus  = lplus.split("=",1)
                            cuts = lplus[1]
                            # print(cutplus, " ++ ", lplus[0])
                            # Check if cut is in file
                            if cutplus in lplus[0]:
                                # Check if run type is already defined in dictionary
                                if typName in cutDict.keys():
                                    if cuts not in cutDict.items():
                                        # If run type already defined, then append dictionary key
                                        if (self.DEBUG):
                                            print("cuts",cuts)
                                        # Grabs parameters from DB (see below)
                                        db_cut = self.search_DB(cuts,runNum)
                                        if (self.DEBUG):
                                            print(typName, " already found!!!!")
                                        cutDict[typName] += ","+db_cut
                                        # print(lplus[0],"++>",cutDict[typName])
                                else:
                                    # If run type not defined, then add key to dictionary
                                    if (self.DEBUG):
                                        print("cuts",cuts)
                                    # Grabs parameters from DB (see below)
                                    db_cut = self.search_DB(cuts,runNum)
                                    cutName = {typName : db_cut}
                                    cutDict.update(cutName)
                                    # print(lplus[0],"++>",cutDict[typName])
                            else:
                                # print("ERROR 6: %s cut does not match %s" % (cutplus,lplus[0]))
                                continue

                    ###################
                    # Subtracted cuts #
                    ###################

                    # Loop over cuts that need to be subtracted
                    for cutminus in minus:
                        if (self.DEBUG):
                            print("- ",cutminus)
                        # Matches run type cuts with the general cuts (e.g pid, track, etc.)
                        if "pid" in cutminus:
                            minusfout = self.REPLAYPATH+"/UTIL_PION/DB/CUTS/general/pid.cuts"
                        elif "track" in cutminus:
                            minusfout = self.REPLAYPATH+"/UTIL_PION/DB/CUTS/general/track.cuts"
                        elif "accept" in cutminus:
                            minusfout = self.REPLAYPATH+"/UTIL_PION/DB/CUTS/general/accept.cuts"
                        elif "coin_time" in cutminus:
                            minusfout = self.REPLAYPATH+"/UTIL_PION/DB/CUTS/general/coin_time.cuts"
                        elif "current" in cutminus:
                            minusfout = self.REPLAYPATH+"/UTIL_PION/DB/CUTS/general/current.cuts"
                        elif "misc" in cutminus:
                            minusfout = self.REPLAYPATH+"/UTIL_PION/DB/CUTS/general/misc.cuts"
                        elif "none" in cutminus:
                            minusfout = "none"
                        else:
                            print("!!!!ERROR!!!!: Subtracted cut %s not defined in /UTIL_PION/DB/CUTS/general/" % cutplus) # ERROR 3
                            print("Cut must be pid, track, accept, coin_time or current")
                            continue
                        # Break down the cut to be removed to find specific leaf to be subtracted from
                        # dictionary
                        minuscut = cutminus.split(".")
                        if len(minuscut) == 3:
                            cutminus = minuscut[1]
                            leafminus = minuscut[2].rstrip()
                        elif minuscut == ['none']:
                            cutminus = "none"
                        else:
                            print("!!!!ERROR!!!!: Invalid syntax for removing cut %s " % (minuscut)) # Error 4
                            continue

                        # Open general cuts file of interest to be removed from dictionary.
                        fminus = open(minusfout)
                            
                        for lminus in fminus:
                            if "#" in lminus:
                                continue
                            else:
                                lminus  = lminus.split("=",1)
                                cuts = lminus[1]
                                # Split cuts to check for the one to be removed.
                                arr_cuts = cuts.split(",")
                                # print(leafminus,": ",cutminus, " -- ", lminus[0])
                                # Check if cut is in file
                                if cutminus in lminus[0]:
                                    for remove in arr_cuts:
                                        # Check which cut matches the one wanted to be removed
                                        if leafminus in remove:
                                            # Grabs parameters from DB (see below)
                                            remove = self.search_DB(remove,runNum)
                                            if (self.DEBUG):
                                                print("Removing... ",remove)
                                            # Replace unwanted cut with blank string
                                            cutDict[typName] = cutDict[typName].replace(remove,"")
                                            # print(lminus[0],"-->",cutDict[typName])
                                else:
                                    # print("ERROR 7: %s cut does not match %s" % (cutminus,lminus[0]))
                                    continue
                        fplus.close()
                        fminus.close()
                if (self.DEBUG):   
                    print("\n\n")
        f.close()
        if (self.DEBUG):
            print(cutDict.keys())
        return cutDict

    # Grabs the cut parameters from the database. In essence this method simply replaces one string
    # with another
    def search_DB(self,cuts,runNum):

        # Split all cuts into a list
        cuts = cuts.split(",")
        db_cuts = []
        for cut in cuts:
            # Find which cut is being called
            if "accept" in cut:
                tmp = cut.split("accept")
                for val in tmp:
                    if "." in val:
                        tmp = val.split(")")[0]
                        tmp = tmp.split(".")[1]
                        fout = self.REPLAYPATH+"/UTIL_PION/DB/PARAM/Acceptance_Parameters.csv"
                        try:
                            data = dict(pd.read_csv(fout))
                        except IOError:
                            print("ERROR 9: %s not found in %s" % (tmp,fout))
                        for i,evt in enumerate(data['Run_Start']):
                            if data['Run_Start'][i] <= np.int64(runNum) <= data['Run_End'][i]:
                                # print("xxxx",tmp,str(data[tmp][i]))
                                cut  = cut.replace("accept."+tmp,str(data[tmp][i]))
                                pass
                            else:
                                # print("!!!!ERROR!!!!: Run %s not found in range %s-%s" % (np.int64(runNum),data['Run_Start'][i],data['Run_End'][i])) # Error 10
                                continue
                    else:
                        continue
                db_cuts.append(cut.rstrip())
            elif "track" in cut:
                tmp = cut.split("track")
                for val in tmp:
                    if "." in val:
                        tmp = val.split(")")[0]
                        tmp = tmp.split(".")[1]
                        fout = self.REPLAYPATH+"/UTIL_PION/DB/PARAM/Tracking_Parameters.csv"
                        try:
                            data = dict(pd.read_csv(fout))
                        except IOError:
                            print("ERROR 9: %s not found in %s" % (tmp,fout))
                        for i,evt in enumerate(data['Run_Start']):
                            if data['Run_Start'][i] <= np.int64(runNum) <= data['Run_End'][i]:
                                cut  = cut.replace("track."+tmp,str(data[tmp][i]))
                                pass
                            else:
                                # print("!!!!ERROR!!!!: Run %s not found in range %s-%s" % (np.int64(runNum),data['Run_Start'][i],data['Run_End'][i])) # Error 10
                                continue
                    else:
                        continue
                db_cuts.append(cut.rstrip())
            elif "CT" in cut:
                tmp = cut.split("CT")
                for val in tmp:
                    if "." in val:
                        tmp = val.split(")")[0]
                        tmp = tmp.split(".")[1]
                        fout = self.REPLAYPATH+"/UTIL_PION/DB/PARAM/Timing_Parameters.csv"
                        try:
                            data = dict(pd.read_csv(fout))
                        except IOError:
                            print("ERROR 9: %s not found in %s" % (tmp,fout))
                        for i,evt in enumerate(data['Run_Start']):
                            if data['Run_Start'][i] <= np.int64(runNum) <= data['Run_End'][i]:
                                cut  = cut.replace("CT."+tmp,str(data[tmp][i]))
                                pass
                            else:
                                # print("!!!!ERROR!!!!: Run %s not found in range %s-%s" % (np.int64(runNum),data['Run_Start'][i],data['Run_End'][i])) # Error 10
                                continue
                    else:
                        continue
                db_cuts.append(cut.rstrip())
            elif "pid" in cut:
                tmp = cut.split("pid")
                for val in tmp:
                    if "." in val:
                        tmp = val.split(")")[0]
                        tmp = tmp.split(".")[1]
                        fout = self.REPLAYPATH+"/UTIL_PION/DB/PARAM/PID_Parameters.csv"
                        try:
                            data = dict(pd.read_csv(fout))
                        except IOError:
                            print("ERROR 9: %s not found in %s" % (tmp,fout))
                        for i,evt in enumerate(data['Run_Start']):
                            if data['Run_Start'][i] <= np.int64(runNum) <= data['Run_End'][i]:
                                cut  = cut.replace("pid."+tmp,str(data[tmp][i]))
                                pass
                            else:
                                # print("!!!!ERROR!!!!: Run %s not found in range %s-%s" % (np.int64(runNum),data['Run_Start'][i],data['Run_End'][i])) # Error 10
                                continue
                    else:
                        continue
                db_cuts.append(cut.rstrip())
            elif "misc" in cut:
                tmp = cut.split("misc")
                for val in tmp:
                    if "." in val:
                        tmp = val.split(")")[0]
                        tmp = tmp.split(".")[1]
                        fout = self.REPLAYPATH+"/UTIL_PION/DB/PARAM/Misc_Parameters.csv"
                        try:
                            data = dict(pd.read_csv(fout))
                        except IOError:
                            print("ERROR 9: %s not found in %s" % (tmp,fout))
                        for i,evt in enumerate(data['Run_Start']):
                            if data['Run_Start'][i] <= np.int64(runNum) <= data['Run_End'][i]:
                                cut  = cut.replace("misc."+tmp,str(data[tmp][i]))
                                pass
                            else:
                                # print("!!!!ERROR!!!!: Run %s not found in range %s-%s" % (np.int64(runNum),data['Run_Start'][i],data['Run_End'][i])) # Error 10
                                continue
                    else:
                        continue
                db_cuts.append(cut.rstrip())                
            # Find which cut is being called. This elif statement is a little different since it only
            # grabs a threshold current value. This is hardcoded for now, eventually need to change.
            elif "current" in cut:
                tmp = cut.split(".")
                tmp = tmp[1].split(")")[0]
                cut  = cut.replace("current."+tmp,"2.5")
                db_cuts.append(cut.rstrip())
                # print("!!!!",cut)
            else:
                # print("ERROR 11: %s not defined" % cut)
                continue
            
        # Rejoins list of cuts to a string separated by commas
        db_cuts  = ','.join(db_cuts)
        return db_cuts

    # Create a working dictionary for cuts by converting string to array of cuts.
    def w_dict(self,cuts):

        inputDict = self.cutDict
        subDict = inputDict[cuts]
        subDict = subDict.split(",")
        cut_arr = [evt for evt in subDict]
        return cut_arr

    # Old version of apply cuts
    def applyCuts(self,leaf,cuts=None):
        
        if cuts:
            tmp = leaf
            applycut = 'tmp['
            i=0
            while i < (len(cuts)-1):
                applycut += 'self.cut("%s") & ' % cuts[i]
                i+=1
            applycut += 'self.cut("%s")]' % cuts[len(cuts)-1]
            tmp = eval(applycut)
        else:
            if (self.DEBUG):
                print('No cuts applied to %s' % leaf)
            tmp = leaf
        
        return tmp

    # New version of applying cuts. The general idea is to apply cuts without sacrificing computation
    # time. Array indexing is much faster than most methods in python. This method formats a string with
    # the cuts required. This string is evaluated and the array index calls the cut() method.See
    # description above for how the analysis script should be formatted. 
    def add_cut(self,arr, cuts):

        arr_cut = arr  
        applycut = "arr_cut["
        inputDict = self.cutDict
        subDict = inputDict[cuts]
        for i,(key,val) in enumerate(subDict.items()):
            if i == len(subDict)-1:
                applycut += 'self.cut("%s","%s")]' % (key,cuts)
            else:
                applycut += 'self.cut("%s","%s") & ' % (key,cuts)
        arr_cut = eval(applycut)        
        return arr_cut

    # The array index that was evaluated in the add_cut() method calls this method. This method then
    # grabs the properly formated dictionary (from class pyDict) and outputs arrays with cuts.
    def cut(self,key,cuts=None):

        if cuts:
            inputDict = self.cutDict
            subDict = inputDict[cuts]
            value = subDict.get(key,"Leaf name not found")
            return value
        # Just for old version for applying cuts (i.e. applyCuts() method)
        else:
            return self.cutDict.get(key,"Leaf name not found")

    # A simple progress bar to use in loops
    def progressBar(self,value, endvalue, bar_length):

        percent = float(value) / endvalue
        arrow = '=' * int(round(percent * bar_length)-1) + '>'
        spaces = ' ' * (bar_length - len(arrow))
        
        sys.stdout.write(" \r[{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
        sys.stdout.flush()

    # Creates nice density plots using matplotlib
    def densityPlot(self,x,y,title,xlabel,ylabel,binx,biny,pyMisc,
                    xmin=None,xmax=None,ymin=None,ymax=None,cuts=None,figure=None,ax=None,layered=True):
        if cuts:
            xcut  = self.applyCuts(x,cuts)
            ycut = self.applyCuts(y,cuts)
        else:
            xcut = x
            ycut = y
        if ax or figure:
            print("")
        else:
            fig, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
        if (xmin or xmax or ymin or ymax):
            # norm=colors.LogNorm() makes colorbar normed and logarithmic
            hist = ax.hist2d(xcut, ycut,bins=(pyMisc.setbin(x,binx,xmin,xmax),pyMisc.setbin(y,biny,ymin,ymax)), norm=colors.LogNorm())
        else:
            # norm=colors.LogNorm() makes colorbar normed and logarithmic
            hist = ax.hist2d(xcut, ycut,bins=(pyMisc.setbin(x,binx),pyMisc.setbin(y,biny)), norm=colors.LogNorm())
        if layered is True :
            plt.colorbar(hist[3], ax=ax, spacing='proportional', label='Number of Events')

        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        inputVal = [x,y]
        
        if (xmin or xmax or ymin or ymax):
            binVal = [pyMisc.setbin(x,binx,xmin,xmax),pyMisc.setbin(y,biny,ymin,ymax)]
        else:
            binVal = [pyMisc.setbin(x,binx),pyMisc.setbin(y,biny)]
        return [binVal, fig]

    # Creates polar plots (useful for kaonlt analysis). Old script, has not been checked in a while.
    def polarPlot(self,theta,r,title,thetalabel,rlabel,bintheta,binr,pyMisc,
                  thetamin=None,thetamax=None,rmin=None,rmax=None,cuts=None,figure=None,ax=None):

        if cuts:
            thetacut  = self.applyCuts(theta,cuts)
            rcut = self.applyCuts(r,cuts)
        else:
            thetacut = theta
            rcut = r
        xy = np.vstack([thetacut, rcut])
        z = stats.gaussian_kde(xy)(xy)
        idx = z.argsort()
        x, y, z = np.array(thetacut)[idx], np.array(rcut)[idx], z[idx]
        if ax or figure:
            # ax = figure.add_subplot(sub,polar=True)
            print("")
        else:
            fig,ax = plt.subplot(111,polar=True)
        if (thetamin or thetamax or rmin or rmax):
            hist = ax.scatter(thetacut, rcut, c=z, edgecolor='', alpha = 0.75)
        else:
            hist = ax.scatter(thetacut, rcut, c=z, edgecolor='', alpha = 0.75)
        ax.grid(True)
        plt.title(title)
        plt.xlabel(thetalabel)
        plt.ylabel(rlabel)
        # plt.colorbar()

        return fig
