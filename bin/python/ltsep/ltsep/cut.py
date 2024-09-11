#! /usr/bin/python

#
# Description: This package will perform many tasks required for l-t separation physics analysis 
# Analysis script required format for applying cuts.
# ================================================================
# Time-stamp: "2024-01-22 15:47:21 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import pandas as pd
import numpy as np
import os

# Import the pathing class of the module. This allows the database chosen based off your current directory. Allows for easier switching between pion and kaon.
from .pathing import SetPath

class SetCuts():
    '''
    SetCuts()

    ----------------------------------------------------------------------------------------------

    This class will set up the cut dictionary as well as 
    apply cuts to arrays.
    '''
    
    def __init__(self, CURRENT_ENV,cutDict=None):
        '''
        __init__(self,CURRENT_ENV,cutDict=None)
                      |           |
                      |           --> cutDict: Sets the dictionary for the class
                      --> CURRENT_ENV: Input current enviroment path

        ----------------------------------------------------------------------------------------------
        
        Constructor of class takes the current enviroment path and an optional dictionary as input
        '''
        self.cutDict = cutDict
        self.REPLAYPATH = SetPath(CURRENT_ENV).getPath("REPLAYPATH")
        self.UTILPATH = SetPath(CURRENT_ENV).getPath("UTILPATH")

    def __str__(self):
        '''
        __str__(self)

        ----------------------------------------------------------------------------------------------

        String representation of class if called as string (eg print(SetCuts))
        '''

        return "{REPLAYPATH : {self.REPLAYPATH}, UTILPATH : {self.UTILPATH}], cutDict : {self.cutDict}}"

    def __repr__(self):
        '''
        __repr__(self)

        ----------------------------------------------------------------------------------------------

        String representation of class if called as is (eg SetCuts)
        '''

        return "SetCuts([{self.REPLAYPATH},{self.UTILPATH}],{self.cutDict})"

    def setbin(self,arr,numbin,xmin=None,xmax=None):
        '''
        setbin(self,arr,numbin,xmin=None,xmax=None)
                    |    |      |         |
                    |    |      |         --> xmax: Upper array set value
                    |    |      --> xmin: Lower array set value
                    |    --> numbin: Number of bins
                    --> arr: Input array
        ----------------------------------------------------------------------------------------------

        A method for defining a bin. This may be called in any matplotlib package plots.
        This will calculate a suitable bin width and use that to equally distribute the bin size
        '''
        
        if (xmin or xmax):
            arr = self.fixBin(arr,xmin,xmax)
        else:
            arr = arr
            
        # Find bin width
        binwidth = (max(abs(arr))-min(abs(arr)))/numbin
        
        # Create array of bin values
        bins = np.arange(min(arr), max(arr) + binwidth, binwidth)

        return bins

    @staticmethod
    def fixBin(arr,low,high):
        '''
        fixBin(self,arr,low,high)
                    |    |   |
                    |    |   --> high: Upper array set value
                    |    --> low: Lower array set value
                    --> arr: Input array to set limits on

        ----------------------------------------------------------------------------------------------

        This method is complimentary to setbin(). This will cut the distribution based off the min 
        and max array values
        '''

        return arr[(arr > low) & (arr < high)]

    def readDict(self,cut,inputDict=None):
        '''
        readDict(self,cut,inputDict=None)
                      |   |
                      |   --> inputDict: Recieves the input Dictionary to be appended with cut names
                      --> cut: Name of cuts to be appended to dictionary keys

        ----------------------------------------------------------------------------------------------

        Reads in cut names and sets them as a key in the general cut dictionary
        '''
        for key,val in self.cutDict.items():
            if key == cut:
                inputDict.update({key : {}})
        return inputDict

    def evalDict(self,cut,eval_xi,inputDict):
        '''
        evalDict(self,cut,eval_xi,inputDict)
                      |   |       |
                      |   |       --> inputDict: Recieves the dictionary created in readDict()
                      |   --> eval_xi: Evaluates the value of booleanDict() to produce cut array
                      --> cut: Cut name being applied to update cut dictionary with properly cut array

        ----------------------------------------------------------------------------------------------

        Updates general cut dictionary with new cut version of arrays as the dictionary value. The 
        eval_xi must be an eval() method to properly apply boolean list comprehension on uncut array.
        '''

        inputDict[cut].update(eval_xi)
        return inputDict

    def importDict(self,inp_cuts,fout,runNum,DEBUG=False):
        '''
        importDict(self,inp_cuts,fout,runNum,DEBUG=False)
                        |        |    |      |
                        |        |    |      --> DEBUG: Debug flag
                        |        |    --> runNum: Run number
                        |        --> fout: file of run type cuts to be applied
                        --> inp_cuts: List of run type cuts to be applied

        ----------------------------------------------------------------------------------------------

        Imports the cut strings and converts them to a dictionary. 
        '''

        # Open run type cuts of interest
        f = open(fout)
        cutDict = {}

        # Matches run type cuts with the general cuts (e.g pid, track, etc.)
        gencutDict = {
            "pid" : self.UTILPATH+"/DB/CUTS/general/pid.cuts",
            "track" : self.UTILPATH+"/DB/CUTS/general/track.cuts",
            "accept" : self.UTILPATH+"/DB/CUTS/general/accept.cuts",
            "coin_time" : self.UTILPATH+"/DB/CUTS/general/coin_time.cuts",
            "CT" : self.UTILPATH+"/DB/CUTS/general/coin_time.cuts",
            "current" : self.UTILPATH+"/DB/CUTS/general/current.cuts",
            "misc" : self.UTILPATH+"/DB/CUTS/general/misc.cuts",
        }

        def genCut(cut_list,add_flag=True):
            '''
            genCut(cut_list,add_flag=True)
                   |        |
                   |        --> add_flag: Flag False if using subtracted cuts
                   --> cut_list: Input list of cuts with ambiguious parameters

            ----------------------------------------------------------------------------------------------

            Function to get the general cuts and calls the search_DB method to get the param values for each cut
            '''
            gencut = []
            # Loop over cut arguments
            for i,val in enumerate(cut_list):
                cutgen = val.split(".")
                if (DEBUG):
                    print("cutgen ", cutgen)
                # Get the general cut name 
                gencut.append(cutgen)
                if cutgen[0].strip() not in gencutDict:
                    print("!!!!ERROR!!!!: Added cut {0} not defined in {1}/DB/CUTS/general/".format(cutgen[0],self.UTILPATH)) # ERROR 2
                    print("Cut must be pid, track, accept, coin_time or current")
            if (DEBUG):
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
                                        if (DEBUG):
                                            print("cuts ",cuts)
                                            print("val ",val)
                                        # Grabs parameters from DB (see below)
                                        db_cut = self.search_DB(cuts,runNum,DEBUG)
                                        if (DEBUG):
                                            print(typName, " already found!!!!")
                                        cutDict[typName] += ","+db_cut
                                else:
                                    # If run type not defined, then add key to dictionary
                                    if (DEBUG):
                                        print("cuts ",cuts)
                                        print("val ",val)
                                    # Grabs parameters from DB (see below)
                                    db_cut = self.search_DB(cuts,runNum,DEBUG)
                                    cutName = {typName : db_cut}
                                    cutDict.update(cutName)
                            else:
                                continue
                        else:
                            # Break down the cut to be removed to find specific leaf to be subtracted from
                            # dictionary
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
                                        remove = self.search_DB(remove,runNum,DEBUG)
                                        if (DEBUG):
                                            print("Removing... ",remove)
                                        # Replace unwanted cut with blank string
                                        cutDict[typName] = cutDict[typName].replace(remove,"")
                f.close()
            return gencut

        def flatten(minus_list):
            '''
            flatten(minus_list)
                    |
                    --> minus_list: Can input any array but used here to flatten minus cut list

            ----------------------------------------------------------------------------------------------

            Flattens multidimensional list
            '''
            flat_list = []
            for e in minus_list:
                if type(e) is list:
                    for i in e:
                        flat_list.append(i)
                else:
                    flat_list.append(e)
            return flat_list

        for ic in inp_cuts:
            if (DEBUG):
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
                        if (DEBUG):
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
                        if (DEBUG):
                            print("+ ", pluscut)
                            print("- ", minuscut)
                        
                        ##############
                        # Added cuts #
                        ##############
                        if (DEBUG):
                            print("Cuts added...")
                        genpluscut = genCut(pluscut)

                        ###################
                        # Subtracted cuts #
                        ###################
                        if (DEBUG):
                            print("Cuts subtracted...")                
                        genminuscut = genCut(minuscut,add_flag=False)
                        break
        f.close()
        if (DEBUG):
            print("\n\n")
            print(cutDict.keys())
            print("\n\n")
        return cutDict

    def search_DB(self, cuts,runNum,DEBUG):
        '''
        search_DB(self,cuts,runNum,DEBUG)
                       |    |      |
                       |    |      --> DEBUG: Debug flag (depends on importDict() flag)
                       |    --> runNum: Run number from importDict()
                       --> cuts: List of added or subtracted cuts from importDict()

        ----------------------------------------------------------------------------------------------

        Grabs the cut parameters from the database. In essence this method simply replaces one string with another
        '''
        # Split all cuts into a list
        cuts = cuts.split(",")
        db_cuts = []
        
        paramDict = {
            "accept" : self.UTILPATH+"/DB/PARAM/Acceptance_Parameters.csv",
            "track" : self.UTILPATH+"/DB/PARAM/Tracking_Parameters.csv",
            "CT" : self.UTILPATH+"/DB/PARAM/Timing_Parameters.csv",
            "pid" : self.UTILPATH+"/DB/PARAM/PID_Parameters.csv",
            "misc" : self.UTILPATH+"/DB/PARAM/Misc_Parameters.csv",
            "current" : self.UTILPATH+"/DB/PARAM/Current_Parameters.csv"
        }

        def grabCutData(paramName,cut):
            '''
            grabCutData(paramName,cut)
                        |         |
                        |         --> cut: General cut who's param file to check
                        --> paramName: Name of parameter column to grab value from

            ----------------------------------------------------------------------------------------------
            
            Grab parameter values from search_DB to replace arbitrary set values in cut strings
            '''
            # Find which cut is being called
            if paramName in cut:
                paramVal = cut.split(paramName)
                for val in paramVal:
                    # Splits string and checks for abs() so that it does not cut string around these curved brackets
                    if "." in val and "abs" not in val:
                        paramVal = val.split(")")[0]
                        paramVal = paramVal.split(".")[1]
                        # Search param dictionary for values based off paramName key
                        fout = paramDict[paramName]
                        try:
                            data = dict(pd.read_csv(fout))
                        except IOError:
                            print("ERROR 9: %s not found in %s" % (paramVal,fout))
                        for i,evt in enumerate(data['Run_Start']):
                            # Check if run number is defined in param file
                            if data['Run_Start'][i] <= np.int64(runNum) <= data['Run_End'][i]:
                                cut  = cut.replace(paramName+"."+paramVal,str(data[paramVal][i]))
                                if (DEBUG):
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

        def has_numbers(inputString):
            '''
            has_numbers(inputString)
                        |
                        --> inputString: String to check if contains a number

            ----------------------------------------------------------------------------------------------
            
            Check if string contains a number. Returns true if number is in string.
            '''
            return any(char.isdigit() for char in inputString)
            
        for cut in cuts:
            # Find which cut is being called
            if "accept." in cut:
                grabCutData("accept",cut)
            elif "track." in cut:
                grabCutData("track",cut)
            elif "CT." in cut:
                grabCutData("CT",cut)
            elif "pid." in cut:
                grabCutData("pid",cut)
            elif "misc." in cut:
                grabCutData("misc",cut)          
            elif "current." in cut:
                grabCutData("current",cut)
            elif has_numbers(cut):
                grabCutData("num",cut)
            else:
                # print("ERROR 11: %s not defined" % cut)
                continue
            
        # Rejoins list of cuts to a string separated by commas
        db_cuts  = ','.join(db_cuts)
        return db_cuts

    def booleanDict(self,cut):
        '''
        booleanDict(self,cut)
                         |
                         --> cut: Input cut to convert to boolean

        ----------------------------------------------------------------------------------------------

        Create a boolean dictionary for cut by converting string to array of pass/no pass cut.
        '''

        inputDict = self.cutDict
        subDict = inputDict[cut]
        subDict = subDict.split(",")
        cut_arr = [evt for evt in subDict]
        cut_arr = list(filter(None,cut_arr)) # Filter out blank string values (ie where cuts were removed)
        return cut_arr

    def apply_cut(self,arr, cut):
        '''
        apply_cut(self,arr, cut)
                     |    |
                     |    --> cut: Run type cut name to impliment to array 
                     --> arr: Input array to be cut

        ----------------------------------------------------------------------------------------------

        Creates the string of cuts to later be evaluated and then converted to boolean array.
        '''

        applycut = "arr["
        inputDict = self.cutDict
        subDict = inputDict[cut]
        for i,(key,val) in enumerate(subDict.items()):
            if i == len(subDict)-1:
                applycut += 'self.cut("%s","%s")]' % (key,cut)
            else:
                applycut += 'self.cut("%s","%s") & ' % (key,cut)
        return applycut

    def add_cut(self, arr, cut):
        '''
        add_cut(self,arr, cut)
                     |    |
                     |    --> cut: Run type cut name to impliment to array 
                     --> arr: Input array to be cut

        ----------------------------------------------------------------------------------------------

        Applies cut. The general idea is to apply cut without sacrificing computation
        time. Array indexing is much faster than most methods in python. This method formats a string with
        the cut required. This string is evaluated and the array index calls the cut() method.See
        description above for how the analysis script should be formatted. 
        '''

        return eval(self.apply_cut(arr, cut))

    def cut(self,key,cut):
        '''
        cut(self,key,cut=None)
                 |   |
                 |   --> cut: Run type cut name to impliment to array (called by add_cut())
                 --> key: Key of cut dictionary to call so it can cut input array called by add_cut()

        ----------------------------------------------------------------------------------------------

        The array index that was evaluated in the add_cut() method calls this method. This method then
        grabs the properly formated dictionary (from class pyDict) and outputs arrays with cut.
        '''

        inputDict = self.cutDict
        subDict = inputDict[cut]
        value = subDict.get(key,"Leaf name not found")
        return value
