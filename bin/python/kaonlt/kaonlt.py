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

cuts = ["<list of cuts>"]

c = klt.pyPlot(None) # See below for pyPlot class definition
readDict = c.read_dict(fout) # read in run type cuts file and makes dictionary

def make_cutDict(cut,inputDict=None):
''
This method calls several methods in kaonlt package. It is required to create properly formated
dictionaries. The evaluation must be in the analysis script because the analysis variables (i.e. the
leaves of interest) are not defined in the kaonlt package. This makes the system more flexible
overall, but a bit more cumbersome in the analysis script. Perhaps one day a better solution will be
implimented.
''
    global c

    c = klt.pyPlot(readDict)
    x = c.w_dict(cut)
    print("\n%s" % cut)
    print(x, "\n")
    
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

for i,c in enumerate(cuts):
    if i == 0:
        cutDict = make_cutDict("%s" % c )
    else:
        cutDict = make_cutDict("%s" % c,cutDict)

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

class pyDict(dict):
    '''
    When calling kaonlt package, you may define a dictionary in the script. This dictionary will contain
    the cuts of interest (defined in a CUTS directory).  These cuts are read in through the read_dict()
    method of the pyPlot() class. The pyDict class is not explicitly called, but rather called implicitly
    by other classes.
    '''
    
    def __init__(self,inputTree):
        self.inputTree = inputTree
        
class pyBranch(pyDict):
    '''
    This class, with its findBranch method, will grab the leaves in a branch using uproot package. This takes the tree as an input.
    '''
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

class pyRoot():
    '''    
    This class is for converting files into root files after the analysis steps
    '''

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

class pyEquation():
    '''            
    This class stores a variety of equations often used in the KaonLT analysis procedure
    '''

    # Define missing mass calculation
    def missmass():
        print("missmass")

class pyPlot(pyDict):
    '''
    This is the most extensive class of the kaonlt package. This class will perform many required tasks
    for doing in depth analysis in python. This class does not require, but will use the pyDict class to
    apply cuts. Set the dictionary to None if no cuts are required.
    '''
    
    def __init__(self, REPLAYPATH,cutDict=None,DEBUG=False):
        self.REPLAYPATH = REPLAYPATH
        self.cutDict = cutDict
        self.DEBUG = DEBUG


    def setbin(self,plot,numbin,xmin=None,xmax=None):
        '''
        A method for defining a bin. This may be called in any matplotlib package plots.
        This will calculate a suitable bin width and use that to equally distribute the bin size
        '''
        
        if (xmin or xmax):
            leaf = self.fixBin(plot,plot,xmin,xmax)
        else:
            leaf = plot
            
        binwidth = (abs(leaf).max()-abs(leaf).min())/numbin
        
        bins = np.arange(min(leaf), max(leaf) + binwidth, binwidth)

        return bins

    def fixBin(self,cut,plot,low,high):
        '''
        This method is complimentary to setbin(). This will cut the distribution based off the min and max array values

        '''
        arrCut = cut
        arrPlot = plot
        arrPlot = arrPlot[(arrCut > low) & (arrCut < high)]

        return arrPlot

    def cut_RF(self,runNum,MaxEvent):
        TimingCutFile = self.REPLAYPATH+'/UTIL_PION/DB/PARAM/Timing_Parameters.csv'
        # rootName = "/lustre19/expphy/volatile/hallc/c-kaonlt/sjdkay/ROOTfiles/Proton_Analysis/Pass3/Proton_coin_replay_production_%s_%s.root" % (self.REPLAYPATH, runNum, MaxEvent)
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

    def read_dict(self,inp_cuts,fout,runNum):
        '''
        This method reads in the CUTS and converts them to a dictionary. 
        '''

        # Open run type cuts of interest
        f = open(fout)
        cutDict = {}

        # Matches run type cuts with the general cuts (e.g pid, track, etc.)
        gencutDict = {
            "pid" : self.REPLAYPATH+"/UTIL_PION/DB/CUTS/general/pid.cuts",
            "track" : self.REPLAYPATH+"/UTIL_PION/DB/CUTS/general/track.cuts",
            "accept" : self.REPLAYPATH+"/UTIL_PION/DB/CUTS/general/accept.cuts",
            "coin_time" : self.REPLAYPATH+"/UTIL_PION/DB/CUTS/general/coin_time.cuts",
            "current" : self.REPLAYPATH+"/UTIL_PION/DB/CUTS/general/current.cuts",
            "misc" : self.REPLAYPATH+"/UTIL_PION/DB/CUTS/general/misc.cuts",
        }

        def genCut(cut_list,add_flag=True):
            '''
            Function to get the general cuts and calls the search_DB method to get the param values for each cut
            '''
            gencut = []
            # Loop over cut arguments
            for i,val in enumerate(cut_list):
                cutgen = val.split(".")
                if (self.DEBUG):
                    print("cutgen ", cutgen)
                # Get the general cut name 
                gencut.append(cutgen)
                if cutgen[0].strip() not in gencutDict:
                    print("!!!!ERROR!!!!: Added cut %s not defined in /UTIL_PION/DB/CUTS/general/" % cutgen[0]) # ERROR 2
                    print("Cut must be pid, track, accept, coin_time or current")
            if (self.DEBUG):
                print("gencuts ", gencut)
            for i,val in enumerate(gencut):
                # Open general cuts file of interest to be added to dictionary
                f = open(gencutDict[val[0]])
                for cutline in f:
                    # Ignore comments
                    if "#" in cutline:
                        continue
                    else:
                        # Redefine cut as 2nd element of list split at = (but only the first instance of =)
                        # This 2nd element are the general cuts
                        cutName = cutline.split("=",1)[0].strip().strip("\n")
                        cuts = cutline.split("=",1)[1].strip().strip("\n")
                        if add_flag:
                            # Check for general cut that was called
                            if val[1] == cutName:
                                # Check if run type is already defined in dictionary
                                if typName in cutDict.keys():
                                    if cuts not in cutDict.items():
                                        # If run type already defined, then append dictionary key
                                        if (self.DEBUG):
                                            print("cuts ",cuts)
                                            print("val ",val)
                                        # Grabs parameters from DB (see below)
                                        db_cut = self.search_DB(cuts,runNum)
                                        if (self.DEBUG):
                                            print(typName, " already found!!!!")
                                        cutDict[typName] += ","+db_cut
                                else:
                                    # If run type not defined, then add key to dictionary
                                    if (self.DEBUG):
                                        print("cuts ",cuts)
                                        print("val ",val)
                                    # Grabs parameters from DB (see below)
                                    db_cut = self.search_DB(cuts,runNum)
                                    cutName = {typName : db_cut}
                                    cutDict.update(cutName)
                            else:
                                continue
                        else:
                            # Break down the cut to be removed to find specific leaf to be subtracted from
                            # dictionary
                            #minuscut = gencut[0]
                            minuscut = val
                            if len(minuscut) == 3:
                                cutminus = minuscut[1]
                                leafminus = minuscut[2].rstrip()
                            elif minuscut == ['none']:
                                cutminus = "none"
                            else:
                                print("!!!!ERROR!!!!: Invalid syntax for removing cut %s " % (minuscut)) # Error 4
                                continue
                            # Split cuts to check for the one to be removed.
                            arr_cuts = cuts.split(",")
                            # Check for general cut that was called
                            if val[1] == cutName:
                                for remove in arr_cuts:
                                    # Check which cut matches the one wanted to be removed
                                    if leafminus in remove:
                                        # Grabs parameters from DB (see below)
                                        remove = self.search_DB(remove,runNum)
                                        if (self.DEBUG):
                                            print("Removing... ",remove)
                                        # Replace unwanted cut with blank string
                                        cutDict[typName] = cutDict[typName].replace(remove,"")
                f.close()
            return gencut

        def flatten(minus_list):
            flat_list = []
            for e in minus_list:
                if type(e) is list:
                    for i in e:
                        flat_list.append(i)
                else:
                    flat_list.append(e)
            return flat_list

        for ic in inp_cuts:
            if (self.DEBUG):
                print("\nInput ", ic)
            f.seek(0)
            for line in f:
                # Ignore comments
                if "#" in line:
                    continue
                else:
                    line = line.split("=",1)
                    # Grab run type cut name
                    typName = line[0].strip()
                    if ic == typName:
                        # Grab run type cuts required, note at this stage the cuts to be removed are bunched
                        # together still
                        pluscut = line[1].split("+")
                        pluscut = [i.strip().strip("\n") for i in pluscut]
                        if (self.DEBUG):
                            print("Type ", typName)
                            print("Cuts ", pluscut)
                        minuscut = [None]*len(pluscut)
                        # Loop over run type cuts being split by +
                        for i,evt in enumerate(pluscut):
                            # Split any cuts to be removed
                            cutminus = evt.split("-")
                            if len(cutminus) > 1:
                                # Define first cut to be added, any other cuts to be added will be done in future
                                # iteration over run type cuts
                                pluscut[i] = cutminus[0].strip()
                                # Ignore first element, since it will always be an added cut
                                minuscut[i] = cutminus[1:]
                        minuscut = flatten([x for x in minuscut if x is not None])
                        if (self.DEBUG):
                            print("+ ", pluscut)
                            print("- ", minuscut)
                        
                        ##############
                        # Added cuts #
                        ##############
                        if (self.DEBUG):
                            print("Cuts added...")
                        genpluscut = genCut(pluscut)

                        ###################
                        # Subtracted cuts #
                        ###################
                        if (self.DEBUG):
                            print("Cuts subtracted...")                
                        genminuscut = genCut(minuscut,add_flag=False)
                        break

        f.close()
        if (self.DEBUG):
            print("\n\n")
            print(cutDict.keys())
            print("\n\n")
        return cutDict

    def search_DB(self,cuts,runNum):
        '''
        Grabs the cut parameters from the database. In essence this method simply replaces one string with another
        '''
        # Split all cuts into a list
        cuts = cuts.split(",")
        db_cuts = []
        
        paramDict = {
            "accept" : self.REPLAYPATH+"/UTIL_PION/DB/PARAM/Acceptance_Parameters.csv",
            "track" : self.REPLAYPATH+"/UTIL_PION/DB/PARAM/Tracking_Parameters.csv",
            "CT" : self.REPLAYPATH+"/UTIL_PION/DB/PARAM/Timing_Parameters.csv",
            "pid" : self.REPLAYPATH+"/UTIL_PION/DB/PARAM/PID_Parameters.csv",
            "misc" : self.REPLAYPATH+"/UTIL_PION/DB/PARAM/Misc_Parameters.csv",
            "current" : self.REPLAYPATH+"/UTIL_PION/DB/PARAM/Current_Parameters.csv"
        }

        def grabCutData(paramName,cut):
            # Find which cut is being called
            if paramName in cut:
                paramVal = cut.split(paramName)
                for val in paramVal:
                    if "." in val and "abs" not in val:
                        paramVal = val.split(")")[0]
                        paramVal = paramVal.split(".")[1]
                        fout = paramDict[paramName]
                        try:
                            data = dict(pd.read_csv(fout))
                        except IOError:
                            print("ERROR 9: %s not found in %s" % (paramVal,fout))
                        for i,evt in enumerate(data['Run_Start']):
                            if data['Run_Start'][i] <= np.int64(runNum) <= data['Run_End'][i]:
                                cut  = cut.replace(paramName+"."+paramVal,str(data[paramVal][i]))
                                if (self.DEBUG):
                                    print("paramVal ",paramVal, "= ",data[paramVal][i])
                                pass
                            else:
                                # print("!!!!ERROR!!!!: Run %s not found in range %s-%s" % (np.int64(runNum),data['Run_Start'][i],data['Run_End'][i])) # Error 10
                                continue
                    else:
                        continue
            if paramName == "num":
                cut = cut
            db_cuts.append(cut.strip())

        # Returns true if number is in string
        def has_numbers(inputString):
            return any(char.isdigit() for char in inputString)
            
        for cut in cuts:
            # Find which cut is being called
            if "accept" in cut:
                grabCutData("accept",cut)
            elif "track" in cut:
                grabCutData("track",cut)
            elif "CT" in cut:
                grabCutData("CT",cut)
            elif "pid" in cut:
                grabCutData("pid",cut)
            elif "misc" in cut:
                grabCutData("misc",cut)          
            elif "current" in cut:
                grabCutData("current",cut)
            elif has_numbers(cut):
                grabCutData("num",cut)
            else:
                # print("ERROR 11: %s not defined" % cut)
                continue
            
        # Rejoins list of cuts to a string separated by commas
        db_cuts  = ','.join(db_cuts)
        return db_cuts

    def w_dict(self,cuts):
        '''
        Create a working dictionary for cuts by converting string to array of cuts.
        '''

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

    def add_cut(self,arr, cuts):
        '''
        New version of applying cuts. The general idea is to apply cuts without sacrificing computation
        time. Array indexing is much faster than most methods in python. This method formats a string with
        the cuts required. This string is evaluated and the array index calls the cut() method.See
        description above for how the analysis script should be formatted. 
        '''

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

    def cut(self,key,cuts=None):
        '''
        The array index that was evaluated in the add_cut() method calls this method. This method then
        grabs the properly formated dictionary (from class pyDict) and outputs arrays with cuts.
        '''

        if cuts:
            inputDict = self.cutDict
            subDict = inputDict[cuts]
            value = subDict.get(key,"Leaf name not found")
            return value
        # Just for old version for applying cuts (i.e. applyCuts() method)
        else:
            return self.cutDict.get(key,"Leaf name not found")

    def progressBar(self,value, endvalue, bar_length):
        '''
        A simple progress bar to use in loops
        '''

        percent = float(value) / endvalue
        arrow = '=' * int(round(percent * bar_length)-1) + '>'
        spaces = ' ' * (bar_length - len(arrow))
        
        sys.stdout.write(" \r[{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
        sys.stdout.flush()

    def densityPlot(self,x,y,title,xlabel,ylabel,binx,biny,pyMisc,
                    xmin=None,xmax=None,ymin=None,ymax=None,cuts=None,figure=None,ax=None,layered=True):
        '''
        Creates nice density plots using matplotlib
        '''
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

    def polarPlot(self,theta,r,title,thetalabel,rlabel,bintheta,binr,pyMisc,
                  thetamin=None,thetamax=None,rmin=None,rmax=None,cuts=None,figure=None,ax=None):
        '''
        Creates polar plots (useful for kaonlt analysis). Old script, has not been checked in a while.
        '''

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
