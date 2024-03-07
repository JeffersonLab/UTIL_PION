#! /usr/bin/python
#
# Description:
# ================================================================
# Time-stamp: "2023-09-11 23:28:12 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
import uproot as up
import sys

from .cut import SetCuts 
from .pathing import SetPath

#########################
# Cython implimentation #
#from .setcut import *  
#########################

class InvalidEntry(Exception):
    '''
    Raise this exception when something goes wrong with the cuts
    '''
    pass

class Root():
    '''    
    Root()

    ----------------------------------------------------------------------------------------------
    ################################################################################################################################################
    \'''
    Define pathing only
    \'''

    # Import package for cuts
    from ltsep import Root

    lt=Root(os.path.realpath(__file__))

    # Add this to all files for more dynamic pathing
    VOLATILEPATH=lt.VOLATILEPATH
    ANALYSISPATH=lt.ANALYSISPATH
    HCANAPATH=lt.HCANAPATH
    REPLAYPATH=lt.REPLAYPATH
    UTILPATH=lt.UTILPATH
    PACKAGEPATH=lt.PACKAGEPATH
    OUTPATH=lt.OUTPATH
    ROOTPATH=lt.ROOTPATH
    REPORTPATH=lt.REPORTPATH
    CUTPATH=lt.CUTPATH
    PARAMPATH=lt.PARAMPATH
    SCRIPTPATH=lt.SCRIPTPATH
    SIMCPATH=lt.SIMCPATH
    LTANAPATH=lt.LTANAPATH
    CACHEPATH=lt.CACHEPATH
    ANATYPE=lt.ANATYPE
    USER=lt.USER
    HOST=lt.HOST
    # Note the OUTPATH is not defined unless RunType argument is given, see below

    # If you wish to explicitly define root branches then do the following...
    import uproot as up
    tree = up.open("<ROOT_FILE_NAME>")["<ROOT_TREE_NAME>"]
    # Convert root leaf to array with uproot
    branch_name  = tree.array("<ROOT_BRANCH_NAME>") # The periods are replaced with underscores

    ################################################################################################################################################
    \'''
    Define pathing with OUTPATH 
    \'''

    # Import package for cuts
    from ltsep import Root

    lt=Root(os.path.realpath(__file__), "<Run Type (HeePCoin, HeePSing_<spec>, SimcCoin, SimcSing, Prod, Plot_<Type>, None)>")

    # Add this to all files for more dynamic pathing
    VOLATILEPATH=lt.VOLATILEPATH
    ANALYSISPATH=lt.ANALYSISPATH
    HCANAPATH=lt.HCANAPATH
    REPLAYPATH=lt.REPLAYPATH
    UTILPATH=lt.UTILPATH
    PACKAGEPATH=lt.PACKAGEPATH
    OUTPATH=lt.OUTPATH
    ROOTPATH=lt.ROOTPATH
    REPORTPATH=lt.REPORTPATH
    CUTPATH=lt.CUTPATH
    PARAMPATH=lt.PARAMPATH
    SCRIPTPATH=lt.SCRIPTPATH
    SIMCPATH=lt.SIMCPATH
    LTANAPATH=lt.LTANAPATH
    CACHEPATH=lt.CACHEPATH
    ANATYPE=lt.ANATYPE
    USER=lt.USER
    HOST=lt.HOST
    OUTPATH=lt.OUTPATH

    # If you wish to explicitly define root branches then do the following...
    import uproot as up
    tree = up.open("<ROOT_FILE_NAME>")["<ROOT_TREE_NAME>"]
    # Convert root leaf to array with uproot
    branch_name  = tree.array("<ROOT_BRANCH_NAME>") # The periods are replaced with underscores

    ################################################################################################################################################
    \'''
    Define pathing with OUTPATH and root branches
    \'''

    # Import package for cuts
    from ltsep import Root

    # Note that now a ROOTPrefix, runNum, and MaxEvent is required
    lt=Root(os.path.realpath(__file__), "<Run Type (HeePCoin, HeePSing_<spec>, SimcCoin, SimcSing, Prod, Plot_<Type>, None)>", ROOTPrefix, runNum, MaxEvent)

    # Add this to all files for more dynamic pathing
    VOLATILEPATH=lt.VOLATILEPATH
    ANALYSISPATH=lt.ANALYSISPATH
    HCANAPATH=lt.HCANAPATH
    REPLAYPATH=lt.REPLAYPATH
    UTILPATH=lt.UTILPATH
    PACKAGEPATH=lt.PACKAGEPATH
    OUTPATH=lt.OUTPATH
    ROOTPATH=lt.ROOTPATH
    REPORTPATH=lt.REPORTPATH
    CUTPATH=lt.CUTPATH
    PARAMPATH=lt.PARAMPATH
    SCRIPTPATH=lt.SCRIPTPATH
    SIMCPATH=lt.SIMCPATH
    LTANAPATH=lt.LTANAPATH
    CACHEPATH=lt.CACHEPATH
    ANATYPE=lt.ANATYPE
    USER=lt.USER
    HOST=lt.HOST
    OUTPATH=lt.OUTPATH

    # This will allow access to a dictionary of root branches depending on the RunType given
    # Note in this example the cut object, c, is only useful for advanced usage. See below for general use.
    # Note the dictionary of cuts as strings, strDict, is a None object as there are no cuts defined.
    proc_root = lt.setup_ana()
    c = proc_root[0] # Cut object
    tree = proc_root[1] # Dictionary of branches
    strDict = proc_root[2] # Dictionary of cuts as strings

    # Call root branches with the dictionary key
    tree['<ROOT_BRANCH_NAME>']

    ################################################################################################################################################
    \'''
    Define pathing with OUTPATH, root branches, and set up cuts
    \'''

    # Import package for cuts
    from ltsep import Root

    # ---> If multple run type files are required then define a new run type file altogether. Do not try to 
    # chain run type files. It can be done, but is computationally wasteful and pointless.
    cut_f = "<path_to_run_type_cut>"

    cuts = ["<runTypeCut1>","<runTypeCut2>",<etc>,...]

    lt=Root(os.path.realpath(__file__), "<Run Type (HeePCoin, HeePSing_<spec>, SimcCoin, SimcSing, Prod, Plot_<Type>, None)>", ROOTPrefix, runNum, MaxEvent, cut_f, cuts)

    # Add this to all files for more dynamic pathing
    VOLATILEPATH=lt.VOLATILEPATH
    ANALYSISPATH=lt.ANALYSISPATH
    HCANAPATH=lt.HCANAPATH
    REPLAYPATH=lt.REPLAYPATH
    UTILPATH=lt.UTILPATH
    PACKAGEPATH=lt.PACKAGEPATH
    OUTPATH=lt.OUTPATH
    ROOTPATH=lt.ROOTPATH
    REPORTPATH=lt.REPORTPATH
    CUTPATH=lt.CUTPATH
    PARAMPATH=lt.PARAMPATH
    SCRIPTPATH=lt.SCRIPTPATH
    SIMCPATH=lt.SIMCPATH
    LTANAPATH=lt.LTANAPATH
    CACHEPATH=lt.CACHEPATH
    ANATYPE=lt.ANATYPE
    USER=lt.USER
    HOST=lt.HOST
    OUTPATH=lt.OUTPATH

    # Arrays are defined in ltsep, no need to redefine.
    # cut_f, cuts are optional flags. If you don't have cuts just leave these blank and the runtype root branches will be accessible, see above.
    # ROOTPrefix is also an optional flag, see above. This means your branches will need to be defined explicitly, see below.
    proc_root = lt.setup_ana()
    c = proc_root[0] # Cut object
    tree = proc_root[1] # Dictionary of branches
    strDict = proc_root[2] # Dictionary of cuts as 

    # Call root branches with the dictionary key
    tree['<ROOT_BRANCH_NAME>']

    # To apply cuts to root branches...
    # c is the cut object used to grab instance of add_cut
    # add_cut() applies the cut, i.e. "<runTypeCut#>", to the branch defined, i.e. tree['<ROOT_BRANCH_NAME>']
    c.add_cut(tree['<ROOT_BRANCH_NAME>'], "<runTypeCut#>")

    ################################################################################################################################################

    # ----> For more info
    from ltsep import Help

    # Some help examples
    Help.info(Root)
    Help.info(SetCuts.importDict)
    Help.path_setup()
    Help.cut_setup()
    Help.searchPathFile(os.path.realpath(__file__))

    ----------------------------------------------------------------------------------------------

    This is the most extensive class of the ltsep package. This class will grab many of the required 
    tasks for doing in depth analysis in python such as define pathing variables and cuts.
    '''

    def __init__(self, CURRENT_ENV, runType="None", ROOTPrefix="", runNum="-1", MaxEvent="-1", cut_f="", cuts=None, DEBUG=False):
        '''
        __init__(self, CURRENT_ENV, ROOTPrefix, runType, runNum, MaxEvent, cut_f, cuts=None, DEBUG=False)
                       |            |           |        |       |         |      |          |
                       |            |           |        |       |         |      |          --> DEBUG: Set true to show debug output
                       |            |           |        |       |         |      --> cuts: Specific cuts in run type cuts file to call
                       |            |           |        |       |         --> cut_f: File of defined run type cuts
                       |            |           |        |       --> MaxEvent: Max number of events replayed
                       |            |           |        --> runNum: Run number
                       |            |           --> runType: Type of run (HeePCoin, HeePSing_<spec>, SimcCoin, SimcSing, Prod, Plot_<Type>, None, etc.)
                       |            --> ROOTPrefix: ROOT prefix as defined by either the Replay script or other analysis scripts
                       --> CURRENT_ENV: Input current enviroment path

        ----------------------------------------------------------------------------------------------
        
        Constructor of class takes the current enviroment path and an optional dictionary as input
        '''
        self.DEBUG = DEBUG
        self.ROOTPrefix = ROOTPrefix
        self.runNum = runNum
        self.MaxEvent = MaxEvent
        self.cuts = cuts
        self.CURRENT_ENV = CURRENT_ENV    
        self.runType = runType

        # Defines dynamic pathing variables
        self.VOLATILEPATH=SetPath(self.CURRENT_ENV).getPath("VOLATILEPATH")
        self.ANALYSISPATH=SetPath(self.CURRENT_ENV).getPath("ANALYSISPATH")
        self.HCANAPATH=SetPath(self.CURRENT_ENV).getPath("HCANAPATH")
        self.REPLAYPATH=SetPath(self.CURRENT_ENV).getPath("REPLAYPATH")
        self.UTILPATH=SetPath(self.CURRENT_ENV).getPath("UTILPATH")
        self.PACKAGEPATH=SetPath(self.CURRENT_ENV).getPath("PACKAGEPATH")
        self.OUTPATH=SetPath(self.CURRENT_ENV).getPath("OUTPATH")
        self.ROOTPATH=SetPath(self.CURRENT_ENV).getPath("ROOTPATH")
        self.REPORTPATH=SetPath(self.CURRENT_ENV).getPath("REPORTPATH")
        self.CUTPATH=SetPath(self.CURRENT_ENV).getPath("CUTPATH")
        self.PARAMPATH=SetPath(self.CURRENT_ENV).getPath("PARAMPATH")
        self.SCRIPTPATH=SetPath(self.CURRENT_ENV).getPath("SCRIPTPATH")
        self.SIMCPATH=SetPath(self.CURRENT_ENV).getPath("SIMCPATH")
        self.LTANAPATH=SetPath(self.CURRENT_ENV).getPath("LTANAPATH")
        self.CACHEPATH=SetPath(self.CURRENT_ENV).getPath("CACHEPATH")
        self.ANATYPE=SetPath(self.CURRENT_ENV).getPath("ANATYPE")
        self.USER=SetPath(self.CURRENT_ENV).getPath("USER")
        self.HOST=SetPath(self.CURRENT_ENV).getPath("HOST",self.DEBUG)

        ################################################################################################################################################
        '''
        Defines Output pathing and cut location
        '''

        self.cut_f = self.UTILPATH+cut_f

        # Add more path setting as needed in a similar manner
        if "HeeP" in self.runType:
            self.OUTPATH = "%s/OUTPUT/Analysis/HeeP" % self.UTILPATH      # Output folder location
        elif "Simc" in self.runType:
            self.OUTPATH = "%s/OUTPUT/Analysis/HeeP" % self.LTANAPATH      # Output folder location
        elif "Prod" in self.runType:
            self.OUTPATH = "%s/OUTPUT/Analysis/%sLT" % (self.UTILPATH,self.ANATYPE)      # Output folder location
        elif "HGCer" in self.runType:
            self.OUTPATH = "%s/OUTPUT/Analysis/%sLT" % (self.UTILPATH,self.ANATYPE)      # Output folder location            
        elif "Hodo" in self.runType:
            self.OUTPATH = "%s/OUTPUT/Calib/Hodo" % self.UTILPATH      # Output folder location
        else:
            self.OUTPATH = "%s/OUTPUT/Analysis/%s" % (self.UTILPATH, self.runType)      # Output folder location
        self.CUTPATH = "%s/DB/CUTS" % self.UTILPATH

        ################################################################################################################################################
        '''
        Check that root/output paths and files exist for use
        '''

        if self.ROOTPrefix is not "":
            if "Plot" in self.runType:
                # Construct the name of the rootfile based upon the info we provided
                if "Prod" in self.runType:
                    self.rootName = "%s/OUTPUT/Analysis/%sLT/%s_%s_%s.root" % (self.UTILPATH, self.ANATYPE, self.runNum, self.MaxEvent, self.ROOTPrefix,)     # Input file location and variables taking
                elif "HGCer" in self.runType:
                    self.rootName = "%s/OUTPUT/Analysis/%sLT/%s_%s_%s.root" % (self.UTILPATH, self.ANATYPE, self.runNum, self.MaxEvent, self.ROOTPrefix,)     # Input file location and variables taking
                elif "HeeP" in self.runType:
                    self.rootName = "%s/OUTPUT/Analysis/HeeP/%s_%s_%s.root" % (self.UTILPATH, self.runNum, self.MaxEvent, self.ROOTPrefix,)     # Input file location and variables taking
                elif "Simc" in self.runType:
                    self.rootName = "%s/OUTPUT/Analysis/HeeP/%s_%s_%s.root" % (self.UTILPATH, self.runNum, self.MaxEvent, self.ROOTPrefix,)     # Input file location and variables taking
                else:
                    self.rootName = "%s/OUTPUT/Analysis/%s/%s_%s_%s.root" % (self.UTILPATH, self.runType, self.runNum, self.MaxEvent, self.ROOTPrefix)     # Input file location and variables taking
                print ("Attempting to process %s" %(self.rootName))
                SetPath(self.CURRENT_ENV).checkDir(self.OUTPATH)
                SetPath(self.CURRENT_ENV).checkFile(self.rootName)
                print("Output path checks out, outputting to %s" % (self.OUTPATH))
            else:
                # Construct the name of the rootfile based upon the info we provided
                if "Prod" in self.runType:
                    self.rootName = "%s/ROOTfiles/Analysis/%sLT/%s_%s_%s.root" % (self.UTILPATH, self.ANATYPE, self.ROOTPrefix, self.runNum, self.MaxEvent)     # Input file location and variables taking
                elif "HGCer" in self.runType:
                    self.rootName = "%s/ROOTfiles/Analysis/%sLT/%s_%s_%s.root" % (self.UTILPATH, self.ANATYPE, self.ROOTPrefix, self.runNum, self.MaxEvent)     # Input file location and variables taking
                elif "Hodo" in self.runType:
                    self.rootName = "%s/ROOTfiles/Analysis/%sLT/%s_%s_%s.root" % (self.UTILPATH, self.ANATYPE, self.ROOTPrefix, self.runNum, self.MaxEvent)     # Input file location and variables taking
                elif "HeeP" in self.runType:
                    self.rootName = "%s/ROOTfiles/Analysis/HeeP/%s_%s_%s.root" % (self.UTILPATH, self.ROOTPrefix, self.runNum, self.MaxEvent)     # Input file location and variables taking
                elif "Simc" in self.runType:
                    self.rootName = "%s/ROOTfiles/Analysis/HeeP/%s_%s_%s.root" % (self.UTILPATH, self.ROOTPrefix, self.runNum, self.MaxEvent)     # Input file location and variables taking
                else:
                    self.rootName = "%s/ROOTfiles/Analysis/%s/%s_%s_%s.root" % (self.UTILPATH, self.runType, self.ROOTPrefix, self.runNum, self.MaxEvent)     # Input file location and variables taking
                print ("Attempting to process %s" %(self.rootName))
                SetPath(self.CURRENT_ENV).checkDir(self.OUTPATH)
                SetPath(self.CURRENT_ENV).checkFile(self.rootName)
                print("Output path checks out, outputting to %s" % (self.OUTPATH))

        ################################################################################################################################################


    def __str__(self):
        '''
        __str__(self)

        ----------------------------------------------------------------------------------------------

        String representation of class if called as string (eg print(SetCuts))
        '''

        return "{REPLAYPATH : {self.REPLAYPATH}, UTILPATH : {self.UTILPATH}}"

    def __repr__(self):
        '''
        __repr__(self)

        ----------------------------------------------------------------------------------------------

        String representation of class if called as is (eg SetCuts)
        '''

        return "Root([{self.REPLAYPATH},{self.UTILPATH}])"  

    def setup_ana(self):
        '''
        This method brings all the data together and makes it accessible to the script. It calls the other 
        methods to define cuts as well as grabs the dictionary of root branches.
        '''

        # Make cut dictionary and convert to boolean list for cut application
        make_cutDict = self.make_cutDict()
        bool_cuts = make_cutDict[0]

        # Get dictionary of branch names
        treeDict = make_cutDict[1]

        # Get dictionary of cut names and values as strings
        strDict = make_cutDict[2]

        return [bool_cuts,treeDict,strDict]

    def check_runType(self):
        '''
        Creates a list of the root branches for a specific run type.
        '''
        
        def_f = "%s/DB/BRANCH_DEF/%sLT/%s" % (self.UTILPATH,self.ANATYPE,self.runType)

        with open(def_f, 'r') as f:
            def_data = f.read().splitlines()
        return def_data

    def make_cutDict(self):
        '''
        This method calls several methods in ltsep package. It is required to create properly formated
        dictionaries. This will define the root branches based off the run type then define the cut object
        which contains the dictionary of cut boolean lists. 
        '''

        # Grab tree from root file
        e_tree = up.open(self.rootName)["T"]

        # Initiate dictionary of root branches
        treeDict = {}

        # 1) Loops over the root branches of a specific run type (defined in UTILPATH/DB/BRANCH_DEF/<RunTypeFile>)
        # 2) Grabs the branch from the root tree (defined above) and defines as array
        # 3) Adds branch to dictionary 
        for branch in self.check_runType():

            # HMS info
            if branch == "H_dc_InsideDipoleExit":
                H_dc_InsideDipoleExit = e_tree.array("H.dc.InsideDipoleExit")    
                treeDict.update({"H_dc_InsideDipoleExit" : H_dc_InsideDipoleExit})
            if branch == "H_hod_goodscinhit":
                H_hod_goodscinhit = e_tree.array("H.hod.goodscinhit")
                treeDict.update({"H_hod_goodscinhit" : H_hod_goodscinhit})
            if branch == "H_hod_goodstarttime":
                H_hod_goodstarttime = e_tree.array("H.hod.goodstarttime")        
                treeDict.update({"H_hod_goodstarttime" : H_hod_goodstarttime})
            if branch == "H_gtr_beta":
                H_gtr_beta = e_tree.array("H.gtr.beta") # Beta is velocity of particle between pairs of hodoscopes
                treeDict.update({"H_gtr_beta" : H_gtr_beta})
            if branch == "H_dc_x_fp":
                H_dc_x_fp = e_tree.array("H.dc.x_fp")
                treeDict.update({"H_dc_x_fp" : H_dc_x_fp})
            if branch == "H_dc_y_fp":
                H_dc_y_fp = e_tree.array("H.dc.y_fp")
                treeDict.update({"H_dc_y_fp" : H_dc_y_fp})
            if branch == "H_dc_xp_fp":
                H_dc_xp_fp = e_tree.array("H.dc.xp_fp") 
                treeDict.update({"H_dc_xp_fp" : H_dc_xp_fp})
            if branch == "H_dc_yp_fp":
                H_dc_yp_fp = e_tree.array("H.dc.yp_fp") 
                treeDict.update({"H_dc_yp_fp" : H_dc_yp_fp})
            if branch == "H_gtr_xp":
                H_gtr_xp = e_tree.array("H.gtr.th")  # xpfp -> Theta
                treeDict.update({"H_gtr_xp" : H_gtr_xp})
            if branch == "H_gtr_yp":
                H_gtr_yp = e_tree.array("H.gtr.ph")  # ypfp -> Phi
                treeDict.update({"H_gtr_yp" : H_gtr_yp})
            if branch == "H_gtr_dp":
                H_gtr_dp = e_tree.array("H.gtr.dp")  # dp is Delta
                treeDict.update({"H_gtr_dp" : H_gtr_dp})
            if branch == "H_gtr_p":
                H_gtr_p = e_tree.array("H.gtr.p")
                treeDict.update({"H_gtr_p" : H_gtr_p})
            if branch == "H_cal_etotnorm":
                H_cal_etotnorm = e_tree.array("H.cal.etotnorm")  
                treeDict.update({"H_cal_etotnorm" : H_cal_etotnorm})
            if branch == "H_cal_etottracknorm":
                H_cal_etottracknorm = e_tree.array("H.cal.etottracknorm")        
                treeDict.update({"H_cal_etottracknorm" : H_cal_etottracknorm})
            if branch == "H_cer_npeSum":
                H_cer_npeSum = e_tree.array("H.cer.npeSum")      
                treeDict.update({"H_cer_npeSum" : H_cer_npeSum})
            if branch == "H_W":
                H_W = e_tree.array("H.kin.primary.W")
                treeDict.update({"H_W" : H_W})
            if branch == "H_cal_etotnorm":
                H_cal_etotnorm = e_tree.array("H.cal.etotnorm")
                treeDict.update({"H_cal_etotnorm" : H_cal_etotnorm})
            if branch == "H_cer_npeSum":
                H_cer_npeSum = e_tree.array("H.cer.npeSum")
                treeDict.update({"H_cer_npeSum" : H_cer_npeSum})
            if branch == "H_gtr_dp":
                H_gtr_dp = e_tree.array("H.gtr.dp")
                treeDict.update({"H_gtr_dp" : H_gtr_dp})
            if branch == "H_tr_tg_th":
                H_tr_tg_th = e_tree.array("H.gtr.th")
                treeDict.update({"H_tr_tg_th" : H_tr_tg_th})
            if branch == "H_tr_tg_ph":
                H_tr_tg_ph = e_tree.array("H.gtr.ph")
                treeDict.update({"H_tr_tg_ph" : H_tr_tg_ph})
            if branch == "H_gtr_beta":
                H_gtr_beta = e_tree.array("H.gtr.beta")
                treeDict.update({"H_gtr_beta" : H_gtr_beta})
            if branch == "H_tr_chi2":
                H_tr_chi2 = e_tree.array("H.tr.chi2")
                treeDict.update({"H_tr_chi2" : H_tr_chi2})
            if branch == "H_tr_ndof":
                H_tr_ndof = e_tree.array("H.tr.ndof")
                treeDict.update({"H_tr_ndof" : H_tr_ndof})
            if branch == "H_hod_goodscinhit":
                H_hod_goodscinhit = e_tree.array("H.hod.goodscinhit")
                treeDict.update({"H_hod_goodscinhit" : H_hod_goodscinhit})
            if branch == "H_hod_betanotrack":
                H_hod_betanotrack = e_tree.array("H.hod.betanotrack")
                treeDict.update({"H_hod_betanotrack" : H_hod_betanotrack})
            if branch == "H_hod_goodstarttime":
                H_hod_goodstarttime = e_tree.array("H.hod.goodstarttime")
                treeDict.update({"H_hod_goodstarttime" : H_hod_goodstarttime})
            if branch == "H_dc_ntrack":
                H_dc_ntrack = e_tree.array("H.dc.ntrack")
                treeDict.update({"H_dc_ntrack" : H_dc_ntrack})
            if branch == "H_dc_1x1_nhit":
                H_dc_1x1_nhit = e_tree.array("H.dc.1x1.nhit")
                treeDict.update({"H_dc_1x1_nhit" : H_dc_1x1_nhit})
            if branch == "H_dc_1u2_nhit":
                H_dc_1u2_nhit = e_tree.array("H.dc.1u2.nhit")
                treeDict.update({"H_dc_1u2_nhit" : H_dc_1u2_nhit})
            if branch == "H_dc_1u1_nhit":
                H_dc_1u1_nhit = e_tree.array("H.dc.1u1.nhit")
                treeDict.update({"H_dc_1u1_nhit" : H_dc_1u1_nhit})
            if branch == "H_dc_1v1_nhit":
                H_dc_1v1_nhit = e_tree.array("H.dc.1v1.nhit")
                treeDict.update({"H_dc_1v1_nhit" : H_dc_1v1_nhit})
            if branch == "H_dc_1x2_nhit":
                H_dc_1x2_nhit = e_tree.array("H.dc.1x2.nhit")
                treeDict.update({"H_dc_1x2_nhit" : H_dc_1x2_nhit})
            if branch == "H_dc_1v2_nhit":
                H_dc_1v2_nhit = e_tree.array("H.dc.1v2.nhit")
                treeDict.update({"H_dc_1v2_nhit" : H_dc_1v2_nhit})
            if branch == "H_dc_2x1_nhit":
                H_dc_2x1_nhit = e_tree.array("H.dc.2x1.nhit")
                treeDict.update({"H_dc_2x1_nhit" : H_dc_2x1_nhit})
            if branch == "H_dc_2u2_nhit":
                H_dc_2u2_nhit = e_tree.array("H.dc.2u2.nhit")
                treeDict.update({"H_dc_2u2_nhit" : H_dc_2u2_nhit})
            if branch == "H_dc_2u1_nhit":
                H_dc_2u1_nhit = e_tree.array("H.dc.2u1.nhit")
                treeDict.update({"H_dc_2u1_nhit" : H_dc_2u1_nhit})
            if branch == "H_dc_2v1_nhit":
                H_dc_2v1_nhit = e_tree.array("H.dc.2v1.nhit")
                treeDict.update({"H_dc_2v1_nhit" : H_dc_2v1_nhit})    
            if branch == "H_dc_2x2_nhit":
                H_dc_2x2_nhit = e_tree.array("H.dc.2x2.nhit")
                treeDict.update({"H_dc_2x2_nhit" : H_dc_2x2_nhit})
            if branch == "H_dc_2v2_nhit":
                H_dc_2v2_nhit = e_tree.array("H.dc.2v2.nhit")
                treeDict.update({"H_dc_2v2_nhit" : H_dc_2v2_nhit})

            # SHMS info
            if branch == "P_cal_fly_earray":
                P_cal_fly_earray = e_tree.array("P.cal.fly.earray")
                treeDict.update({"P_cal_fly_earray" : P_cal_fly_earray})
            if branch == "P_cal_pr_eplane":
                P_cal_pr_eplane = e_tree.array("P.cal.pr.eplane")
                treeDict.update({"P_cal_pr_eplane" : P_cal_pr_eplane})
            if branch == "P_cal_etotnorm":
                P_cal_etotnorm = e_tree.array("P.cal.etotnorm")
                treeDict.update({"P_cal_etotnorm" : P_cal_etotnorm})
            if branch == "P_aero_npeSum":
                P_aero_npeSum = e_tree.array("P.aero.npeSum")
                treeDict.update({"P_aero_npeSum" : P_aero_npeSum})
            if branch == "P_hgcer_npeSum":
                P_hgcer_npeSum = e_tree.array("P.hgcer.npeSum")
                treeDict.update({"P_hgcer_npeSum" : P_hgcer_npeSum})
            if branch == "P_hgcer_xAtCer":
                P_hgcer_xAtCer = e_tree.array("P.hgcer.xAtCer")
                treeDict.update({"P_hgcer_xAtCer" : P_hgcer_xAtCer})
            if branch == "P_hgcer_yAtCer":
                P_hgcer_yAtCer = e_tree.array("P.hgcer.yAtCer")
                treeDict.update({"P_hgcer_yAtCer" : P_hgcer_yAtCer})
            if branch == "P_aero_xAtCer":
                P_aero_xAtCer = e_tree.array("P.aero.xAtAero")
                treeDict.update({"P_aero_xAtCer" : P_aero_xAtCer})
            if branch == "P_aero_yAtCer":
                P_aero_yAtCer = e_tree.array("P.aero.yAtAero")
                treeDict.update({"P_aero_yAtCer" : P_aero_yAtCer})
            if branch == "P_dc_InsideDipoleExit":
                P_dc_InsideDipoleExit = e_tree.array("P.dc.InsideDipoleExit")    
                treeDict.update({"P_dc_InsideDipoleExit" : P_dc_InsideDipoleExit})
            if branch == "P_hod_goodscinhit":    
                P_hod_goodscinhit = e_tree.array("P.hod.goodscinhit")
                treeDict.update({"P_hod_goodscinhit" : P_hod_goodscinhit})
            if branch == "P_hod_goodstarttime":
                P_hod_goodstarttime = e_tree.array("P.hod.goodstarttime")        
                treeDict.update({"P_hod_goodstarttime" : P_hod_goodstarttime})
            if branch == "P_gtr_beta":
                P_gtr_beta = e_tree.array("P.gtr.beta") # Beta is velocity of particle between pairs of hodoscopes
                treeDict.update({"P_gtr_beta" : P_gtr_beta})
            if branch == "P_gtr_x":
                P_gtr_x = e_tree.array("P.gtr.x")
                treeDict.update({"P_gtr_x" : P_gtr_x})
            if branch == "P_gtr_y":
                P_gtr_y = e_tree.array("P.gtr.y")
                treeDict.update({"P_gtr_y" : P_gtr_y})       
            if branch == "P_dc_x_fp":
                P_dc_x_fp = e_tree.array("P.dc.x_fp")
                treeDict.update({"P_dc_x_fp" : P_dc_x_fp})
            if branch == "P_dc_y_fp":
                P_dc_y_fp = e_tree.array("P.dc.y_fp")
                treeDict.update({"P_dc_y_fp" : P_dc_y_fp})
            if branch == "P_dc_xp_fp":
                P_dc_xp_fp = e_tree.array("P.dc.xp_fp") 
                treeDict.update({"P_dc_xp_fp" : P_dc_xp_fp})
            if branch == "P_dc_yp_fp":
                P_dc_yp_fp = e_tree.array("P.dc.yp_fp") 
                treeDict.update({"P_dc_yp_fp" : P_dc_yp_fp})
            if branch == "P_gtr_xp":
                P_gtr_xp = e_tree.array("P.gtr.th")  # xpfp -> Theta
                treeDict.update({"P_gtr_xp" : P_gtr_xp})
            if branch == "P_gtr_yp":
                P_gtr_yp = e_tree.array("P.gtr.ph")  # ypfp -> Phi
                treeDict.update({"P_gtr_yp" : P_gtr_yp})
            if branch == "P_gtr_p":
                P_gtr_p = e_tree.array("P.gtr.p")
                treeDict.update({"P_gtr_p" : P_gtr_p})
            if branch == "P_gtr_dp":
                P_gtr_dp = e_tree.array("P.gtr.dp")  # dp is Delta
                treeDict.update({"P_gtr_dp" : P_gtr_dp})
            if branch == "P_cal_etotnorm":
                P_cal_etotnorm = e_tree.array("P.cal.etotnorm")  
                treeDict.update({"P_cal_etotnorm" : P_cal_etotnorm})
            if branch == "P_cal_etottracknorm":
                P_cal_etottracknorm = e_tree.array("P.cal.etottracknorm")        
                treeDict.update({"P_cal_etottracknorm" : P_cal_etottracknorm})
            if branch == "P_aero_npeSum":
                P_aero_npeSum = e_tree.array("P.aero.npeSum")    
                treeDict.update({"P_aero_npeSum" : P_aero_npeSum})
            if branch == "P_aero_xAtAero":
                P_aero_xAtAero = e_tree.array("P.aero.xAtAero")  
                treeDict.update({"P_aero_xAtAero" : P_aero_xAtAero})
            if branch == "P_aero_yAtAero":
                P_aero_yAtAero = e_tree.array("P.aero.yAtAero")  
                treeDict.update({"P_aero_yAtAero" : P_aero_yAtAero})
            if branch == "P_hgcer_npeSum":
                P_hgcer_npeSum = e_tree.array("P.hgcer.npeSum")  
                treeDict.update({"P_hgcer_npeSum" : P_hgcer_npeSum})
            if branch == "P_hgcer_xAtCer":
                P_hgcer_xAtCer = e_tree.array("P.hgcer.xAtCer")  
                treeDict.update({"P_hgcer_xAtCer" : P_hgcer_xAtCer})
            if branch == "P_hgcer_yAtCer":
                P_hgcer_yAtCer = e_tree.array("P.hgcer.yAtCer")  
                treeDict.update({"P_hgcer_yAtCer" : P_hgcer_yAtCer})
            if branch == "P_cal_etotnorm":
                P_cal_etotnorm = e_tree.array("P.cal.etotnorm")
                treeDict.update({"P_cal_etotnorm" : P_cal_etotnorm})
            if branch == "P_hgcer_npeSum":
                P_hgcer_npeSum = e_tree.array("P.hgcer.npeSum")
                treeDict.update({"P_hgcer_npeSum" : P_hgcer_npeSum})
            if branch == "P_aero_npeSum":
                P_aero_npeSum = e_tree.array("P.aero.npeSum")
                treeDict.update({"P_aero_npeSum" : P_aero_npeSum})
            if branch == "P_gtr_dp":
                P_gtr_dp = e_tree.array("P.gtr.dp")
                treeDict.update({"P_gtr_dp" : P_gtr_dp})
            if branch == "P_gtr_th":
                P_gtr_th = e_tree.array("P.gtr.th")
                treeDict.update({"P_gtr_th" : P_gtr_th})
            if branch == "P_gtr_ph":
                P_gtr_ph = e_tree.array("P.gtr.ph")
                treeDict.update({"P_gtr_ph" : P_gtr_ph})
            if branch == "P_gtr_beta":
                P_gtr_beta = e_tree.array("P.gtr.beta")
                treeDict.update({"P_gtr_beta" : P_gtr_beta})
            if branch == "P_tr_chi2":
                P_tr_chi2 = e_tree.array("P.tr.chi2")
                treeDict.update({"P_tr_chi2" : P_tr_chi2})
            if branch == "P_tr_ndof":
                P_tr_ndof = e_tree.array("P.tr.ndof")
                treeDict.update({"P_tr_ndof" : P_tr_ndof})
            if branch == "P_hod_goodscinhit":
                P_hod_goodscinhit = e_tree.array("P.hod.goodscinhit")
                treeDict.update({"P_hod_goodscinhit" : P_hod_goodscinhit})
            if branch == "P_hod_betanotrack":
                P_hod_betanotrack = e_tree.array("P.hod.betanotrack")
                treeDict.update({"P_hod_betanotrack" : P_hod_betanotrack})
            if branch == "P_hod_goodstarttime":
                P_hod_goodstarttime = e_tree.array("P.hod.goodstarttime")
                treeDict.update({"P_hod_goodstarttime" : P_hod_goodstarttime})
            if branch == "P_dc_ntrack":
                P_dc_ntrack = e_tree.array("P.dc.ntrack")
                treeDict.update({"P_dc_ntrack" : P_dc_ntrack})
            if branch == "P_dc_1x1_nhit":
                P_dc_1x1_nhit = e_tree.array("P.dc.1x1.nhit")
                treeDict.update({"P_dc_1x1_nhit" : P_dc_1x1_nhit})
            if branch == "P_dc_1u2_nhit":
                P_dc_1u2_nhit = e_tree.array("P.dc.1u2.nhit")
                treeDict.update({"P_dc_1u2_nhit" : P_dc_1u2_nhit})
            if branch == "P_dc_1u1_nhit":
                P_dc_1u1_nhit = e_tree.array("P.dc.1u1.nhit")
                treeDict.update({"P_dc_1u1_nhit" : P_dc_1u1_nhit})
            if branch == "P_dc_1v1_nhit":
                P_dc_1v1_nhit = e_tree.array("P.dc.1v1.nhit")
                treeDict.update({"P_dc_1v1_nhit" : P_dc_1v1_nhit})
            if branch == "P_dc_1x2_nhit":
                P_dc_1x2_nhit = e_tree.array("P.dc.1x2.nhit")
                treeDict.update({"P_dc_1x2_nhit" : P_dc_1x2_nhit})
            if branch == "P_dc_1v2_nhit":
                P_dc_1v2_nhit = e_tree.array("P.dc.1v2.nhit")
                treeDict.update({"P_dc_1v2_nhit" : P_dc_1v2_nhit})
            if branch == "P_dc_2x1_nhit":
                P_dc_2x1_nhit = e_tree.array("P.dc.2x1.nhit")
                treeDict.update({"P_dc_2x1_nhit" : P_dc_2x1_nhit})
            if branch == "P_dc_2u2_nhit":
                P_dc_2u2_nhit = e_tree.array("P.dc.2u2.nhit")
                treeDict.update({"P_dc_2u2_nhit" : P_dc_2u2_nhit})
            if branch == "P_dc_2u1_nhit":
                P_dc_2u1_nhit = e_tree.array("P.dc.2u1.nhit")
                treeDict.update({"P_dc_2u1_nhit" : P_dc_2u1_nhit})       
            if branch == "P_dc_2v1_nhit":
                P_dc_2v1_nhit = e_tree.array("P.dc.2v1.nhit")
                treeDict.update({"P_dc_2v1_nhit" : P_dc_2v1_nhit})
            if branch == "P_dc_2x2_nhit":
                P_dc_2x2_nhit = e_tree.array("P.dc.2x2.nhit")
                treeDict.update({"P_dc_2x2_nhit" : P_dc_2x2_nhit})
            if branch == "P_dc_2v2_nhit":
                P_dc_2v2_nhit = e_tree.array("P.dc.2v2.nhit")
                treeDict.update({"P_dc_2v2_nhit" : P_dc_2v2_nhit})

            # Raster
            if branch == "raster_x":
                raster_x = e_tree.array("P.rb.x")
                treeDict.update({"raster_x" : raster_x})
            if branch == "raster_y":
                raster_y = e_tree.array("P.rb.y")
                treeDict.update({"raster_y" : raster_y})
            if branch == "raster_z":
                raster_z = e_tree.array("P.rb.z")
                treeDict.update({"raster_z" : raster_z})

            # BPM target
            if branch == "bpm_tar_x":
                bpm_tar_x = e_tree.array("P.rb.raster.fr_xbpm_tar")
                treeDict.update({"bpm_tar_x" : bpm_tar_x})
            if branch == "bpm_tar_y":
                bpm_tar_y = e_tree.array("P.rb.raster.fr_ybpm_tar")
                treeDict.update({"bpm_tar_y" : bpm_tar_y})                
                
            # Kinematic quantitites
            if branch == "Q2":
                Q2 = e_tree.array("H.kin.primary.Q2")
                treeDict.update({"Q2" : Q2})
            if branch == "W":
                W = e_tree.array("H.kin.primary.W")  
                treeDict.update({"W" : W})
            if branch == "epsilon":
                epsilon = e_tree.array("H.kin.primary.epsilon")  
                treeDict.update({"epsilon" : epsilon})
            if branch == "ph_q":
                ph_q = e_tree.array("P.kin.secondary.ph_xq")     
                treeDict.update({"ph_q" : ph_q})
            if branch == "ph_recoil":
                ph_recoil = e_tree.array("P.kin.secondary.ph_bq")     
                treeDict.update({"ph_recoil" : ph_recoil})
            if branch == "th_q":
                th_q = e_tree.array("P.kin.secondary.th_xq")     
                treeDict.update({"th_q" : th_q})
            if branch == "th_recoil":
                th_recoil = e_tree.array("P.kin.secondary.th_bq")     
                treeDict.update({"th_recoil" : th_recoil})
            if branch == "emiss":
                emiss = e_tree.array("P.kin.secondary.emiss")    
                treeDict.update({"emiss" : emiss})
            if branch == "MMpi":
                MMpi = e_tree.array("P.kin.secondary.MMpi")      
                treeDict.update({"MMpi" : MMpi})
            if branch == "MMK":
                MMK = e_tree.array("P.kin.secondary.MMK")        
                treeDict.update({"MMK" : MMK})
            if branch == "MMp":
                MMp = e_tree.array("P.kin.secondary.MMp")        
                treeDict.update({"MMp" : MMp})
            if branch == "MandelT":
                MandelT = e_tree.array("P.kin.secondary.MandelT")
                treeDict.update({"MandelT" : MandelT})
            if branch == "MandelU":
                MandelU = e_tree.array("P.kin.secondary.MandelU")   
                treeDict.update({"MandelU" : MandelU})
            if branch == "pmiss":
                pmiss = e_tree.array("P.kin.secondary.pmiss")    
                treeDict.update({"pmiss" : pmiss})
            if branch == "pmiss_x":
                pmiss_x = e_tree.array("P.kin.secondary.pmiss_x")
                treeDict.update({"pmiss_x" : pmiss_x})
            if branch == "pmiss_y":
                pmiss_y = e_tree.array("P.kin.secondary.pmiss_y")
                treeDict.update({"pmiss_y" : pmiss_y})
            if branch == "pmiss_z":
                pmiss_z = e_tree.array("P.kin.secondary.pmiss_z")
                treeDict.update({"pmiss_z" : pmiss_z})
            if branch == "Erecoil":
                Erecoil = e_tree.array("P.kin.secondary.Erecoil")
                treeDict.update({"Erecoil" : Erecoil})
            if branch == "emiss_nuc":
                emiss_nuc = e_tree.array("P.kin.secondary.emiss_nuc")
                treeDict.update({"emiss_nuc" : emiss_nuc})
            if branch == "Mrecoil":
                Mrecoil = e_tree.array("P.kin.secondary.Mrecoil")
                treeDict.update({"Mrecoil" : Mrecoil})                

            # Current
            if branch == "H_bcm_bcm1_AvgCurrent":
                H_bcm_bcm1_AvgCurrent = e_tree.array("H.bcm.bcm1.AvgCurrent")
                treeDict.update({"H_bcm_bcm1_AvgCurrent" : H_bcm_bcm1_AvgCurrent})
            if branch == "H_bcm_bcm2_AvgCurrent":
                H_bcm_bcm2_AvgCurrent = e_tree.array("H.bcm.bcm2.AvgCurrent")
                treeDict.update({"H_bcm_bcm2_AvgCurrent" : H_bcm_bcm2_AvgCurrent})
            if branch == "H_bcm_bcm4a_AvgCurrent":
                H_bcm_bcm4a_AvgCurrent = e_tree.array("H.bcm.bcm4a.AvgCurrent")
                treeDict.update({"H_bcm_bcm4a_AvgCurrent" : H_bcm_bcm4a_AvgCurrent})
            if branch == "H_bcm_bcm4b_AvgCurrent":
                H_bcm_bcm4b_AvgCurrent = e_tree.array("H.bcm.bcm4b.AvgCurrent")
                treeDict.update({"H_bcm_bcm4b_AvgCurrent" : H_bcm_bcm4b_AvgCurrent})
            if branch == "H_bcm_bcm4c_AvgCurrent":
                H_bcm_bcm4c_AvgCurrent = e_tree.array("H.bcm.bcm4c.AvgCurrent")
                treeDict.update({"H_bcm_bcm4c_AvgCurrent" : H_bcm_bcm4c_AvgCurrent})

            # Timing info
            if branch == "CTime_eKCoinTime_ROC1":
                CTime_eKCoinTime_ROC1 = e_tree.array("CTime.eKCoinTime_ROC1")
                treeDict.update({"CTime_eKCoinTime_ROC1" : CTime_eKCoinTime_ROC1})
            if branch == "CTime_ePiCoinTime_ROC1":
                CTime_ePiCoinTime_ROC1 = e_tree.array("CTime.ePiCoinTime_ROC1")
                treeDict.update({"CTime_ePiCoinTime_ROC1" : CTime_ePiCoinTime_ROC1})
            if branch == "CTime_epCoinTime_ROC1":
                CTime_epCoinTime_ROC1 = e_tree.array("CTime.epCoinTime_ROC1")
                treeDict.update({"CTime_epCoinTime_ROC1" : CTime_epCoinTime_ROC1})

            if branch == "P_RF_tdcTime":
                P_RF_tdcTime = e_tree.array("T.coin.pRF_tdcTime")  
                treeDict.update({"P_RF_tdcTime" : P_RF_tdcTime})
            if branch == "P_hod_fpHitsTime":
                P_hod_fpHitsTime = e_tree.array("P.hod.fpHitsTime")
                treeDict.update({"P_hod_fpHitsTime" : P_hod_fpHitsTime})
            if branch == "H_RF_Dist":
                H_RF_Dist = e_tree.array("RFTime.HMS_RFtimeDist")
                treeDict.update({"H_RF_Dist" : H_RF_Dist})
            if branch == "P_RF_Dist":
                P_RF_Dist = e_tree.array("RFTime.SHMS_RFtimeDist")  
                treeDict.update({"P_RF_Dist" : P_RF_Dist})

                
            if branch == "T_coin_pTRIG1_ROC1_tdcTimeRaw":
                T_coin_pTRIG1_ROC1_tdcTimeRaw = e_tree.array("T.coin.pTRIG1_ROC1_tdcTimeRaw")
                treeDict.update({"T_coin_pTRIG1_ROC1_tdcTimeRaw" : T_coin_pTRIG1_ROC1_tdcTimeRaw})
            if branch == "T_coin_pTRIG1_ROC2_tdcTimeRaw":
                T_coin_pTRIG1_ROC2_tdcTimeRaw = e_tree.array("T.coin.pTRIG1_ROC2_tdcTimeRaw")
                treeDict.update({"T_coin_pTRIG1_ROC2_tdcTimeRaw" : T_coin_pTRIG1_ROC2_tdcTimeRaw})
            if branch == "T_coin_pTRIG1_ROC1_tdcTime":
                T_coin_pTRIG1_ROC1_tdcTime = e_tree.array("T.coin.pTRIG1_ROC1_tdcTime")
                treeDict.update({"T_coin_pTRIG1_ROC1_tdcTime" : T_coin_pTRIG1_ROC1_tdcTime})
            if branch == "T_coin_pTRIG1_ROC2_tdcTime":
                T_coin_pTRIG1_ROC2_tdcTime = e_tree.array("T.coin.pTRIG1_ROC2_tdcTime")
                treeDict.update({"T_coin_pTRIG1_ROC2_tdcTime" : T_coin_pTRIG1_ROC2_tdcTime})

            if branch == "T_coin_pTRIG2_ROC1_tdcTimeRaw":
                T_coin_pTRIG2_ROC1_tdcTimeRaw = e_tree.array("T.coin.pTRIG2_ROC1_tdcTimeRaw")
                treeDict.update({"T_coin_pTRIG2_ROC1_tdcTimeRaw" : T_coin_pTRIG2_ROC1_tdcTimeRaw})
            if branch == "T_coin_pTRIG2_ROC2_tdcTimeRaw":
                T_coin_pTRIG2_ROC2_tdcTimeRaw = e_tree.array("T.coin.pTRIG2_ROC2_tdcTimeRaw")
                treeDict.update({"T_coin_pTRIG2_ROC2_tdcTimeRaw" : T_coin_pTRIG2_ROC2_tdcTimeRaw})
            if branch == "T_coin_pTRIG2_ROC1_tdcTime":
                T_coin_pTRIG2_ROC1_tdcTime = e_tree.array("T.coin.pTRIG2_ROC1_tdcTime")
                treeDict.update({"T_coin_pTRIG2_ROC1_tdcTime" : T_coin_pTRIG2_ROC1_tdcTime})
            if branch == "T_coin_pTRIG2_ROC2_tdcTime":
                T_coin_pTRIG2_ROC2_tdcTime = e_tree.array("T.coin.pTRIG2_ROC2_tdcTime")
                treeDict.update({"T_coin_pTRIG2_ROC2_tdcTime" : T_coin_pTRIG2_ROC2_tdcTime})

            if branch == "T_coin_pTRIG3_ROC1_tdcTimeRaw":
                T_coin_pTRIG3_ROC1_tdcTimeRaw = e_tree.array("T.coin.pTRIG3_ROC1_tdcTimeRaw")
                treeDict.update({"T_coin_pTRIG3_ROC1_tdcTimeRaw" : T_coin_pTRIG3_ROC1_tdcTimeRaw})
            if branch == "T_coin_pTRIG3_ROC2_tdcTimeRaw":
                T_coin_pTRIG3_ROC2_tdcTimeRaw = e_tree.array("T.coin.pTRIG3_ROC2_tdcTimeRaw")
                treeDict.update({"T_coin_pTRIG3_ROC2_tdcTimeRaw" : T_coin_pTRIG3_ROC2_tdcTimeRaw})
            if branch == "T_coin_pTRIG3_ROC1_tdcTime":
                T_coin_pTRIG3_ROC1_tdcTime = e_tree.array("T.coin.pTRIG3_ROC1_tdcTime")
                treeDict.update({"T_coin_pTRIG3_ROC1_tdcTime" : T_coin_pTRIG3_ROC1_tdcTime})
            if branch == "T_coin_pTRIG3_ROC2_tdcTime":
                T_coin_pTRIG3_ROC2_tdcTime = e_tree.array("T.coin.pTRIG3_ROC2_tdcTime")
                treeDict.update({"T_coin_pTRIG3_ROC2_tdcTime" : T_coin_pTRIG3_ROC2_tdcTime})

            if branch == "T_coin_pTRIG4_ROC1_tdcTimeRaw":
                T_coin_pTRIG4_ROC1_tdcTimeRaw = e_tree.array("T.coin.pTRIG4_ROC1_tdcTimeRaw")
                treeDict.update({"T_coin_pTRIG4_ROC1_tdcTimeRaw" : T_coin_pTRIG4_ROC1_tdcTimeRaw})
            if branch == "T_coin_pTRIG4_ROC2_tdcTimeRaw":
                T_coin_pTRIG4_ROC2_tdcTimeRaw = e_tree.array("T.coin.pTRIG4_ROC2_tdcTimeRaw")
                treeDict.update({"T_coin_pTRIG4_ROC2_tdcTimeRaw" : T_coin_pTRIG4_ROC2_tdcTimeRaw})
            if branch == "T_coin_pTRIG4_ROC1_tdcTime":
                T_coin_pTRIG4_ROC1_tdcTime = e_tree.array("T.coin.pTRIG4_ROC1_tdcTime")
                treeDict.update({"T_coin_pTRIG4_ROC1_tdcTime" : T_coin_pTRIG4_ROC1_tdcTime})
            if branch == "T_coin_pTRIG4_ROC2_tdcTime":
                T_coin_pTRIG4_ROC2_tdcTime = e_tree.array("T.coin.pTRIG4_ROC2_tdcTime")
                treeDict.update({"T_coin_pTRIG4_ROC2_tdcTime" : T_coin_pTRIG4_ROC2_tdcTime})

            if branch == "T_coin_pTRIG5_ROC1_tdcTimeRaw":
                T_coin_pTRIG5_ROC1_tdcTimeRaw = e_tree.array("T.coin.pTRIG5_ROC1_tdcTimeRaw")
                treeDict.update({"T_coin_pTRIG5_ROC1_tdcTimeRaw" : T_coin_pTRIG5_ROC1_tdcTimeRaw})
            if branch == "T_coin_pTRIG5_ROC2_tdcTimeRaw":
                T_coin_pTRIG5_ROC2_tdcTimeRaw = e_tree.array("T.coin.pTRIG5_ROC2_tdcTimeRaw")
                treeDict.update({"T_coin_pTRIG5_ROC2_tdcTimeRaw" : T_coin_pTRIG5_ROC2_tdcTimeRaw})
            if branch == "T_coin_pTRIG5_ROC1_tdcTime":
                T_coin_pTRIG5_ROC1_tdcTime = e_tree.array("T.coin.pTRIG5_ROC1_tdcTime")
                treeDict.update({"T_coin_pTRIG5_ROC1_tdcTime" : T_coin_pTRIG5_ROC1_tdcTime})
            if branch == "T_coin_pTRIG5_ROC2_tdcTime":
                T_coin_pTRIG5_ROC2_tdcTime = e_tree.array("T.coin.pTRIG5_ROC2_tdcTime")
                treeDict.update({"T_coin_pTRIG5_ROC2_tdcTime" : T_coin_pTRIG5_ROC2_tdcTime})

            if branch == "T_coin_pTRIG6_ROC1_tdcTimeRaw":
                T_coin_pTRIG6_ROC1_tdcTimeRaw = e_tree.array("T.coin.pTRIG6_ROC1_tdcTimeRaw")
                treeDict.update({"T_coin_pTRIG6_ROC1_tdcTimeRaw" : T_coin_pTRIG6_ROC1_tdcTimeRaw})
            if branch == "T_coin_pTRIG6_ROC2_tdcTimeRaw":
                T_coin_pTRIG6_ROC2_tdcTimeRaw = e_tree.array("T.coin.pTRIG6_ROC2_tdcTimeRaw")
                treeDict.update({"T_coin_pTRIG6_ROC2_tdcTimeRaw" : T_coin_pTRIG6_ROC2_tdcTimeRaw})
            if branch == "T_coin_pTRIG6_ROC1_tdcTime":
                T_coin_pTRIG6_ROC1_tdcTime = e_tree.array("T.coin.pTRIG6_ROC1_tdcTime")
                treeDict.update({"T_coin_pTRIG6_ROC1_tdcTime" : T_coin_pTRIG6_ROC1_tdcTime})
            if branch == "T_coin_pTRIG6_ROC2_tdcTime":
                T_coin_pTRIG6_ROC2_tdcTime = e_tree.array("T.coin.pTRIG6_ROC2_tdcTime")
                treeDict.update({"T_coin_pTRIG6_ROC2_tdcTime" : T_coin_pTRIG6_ROC2_tdcTime})

            if branch == "T_coin_pFADC_TREF_ROC2_adcPed":
                T_coin_pFADC_TREF_ROC2_adcPed = e_tree.array("T.coin.pFADC_TREF_ROC2_adcPed")
                treeDict.update({"T_coin_pFADC_TREF_ROC2_adcPed" : T_coin_pFADC_TREF_ROC2_adcPed})
            if branch == "T_coin_hFADC_TREF_ROC1_adcPed":
                T_coin_hFADC_TREF_ROC1_adcPed = e_tree.array("T.coin.hFADC_TREF_ROC1_adcPed")
                treeDict.update({"T_coin_hFADC_TREF_ROC1_adcPed" : T_coin_hFADC_TREF_ROC1_adcPed})
            if branch == "T_coin_pFADC_TREF_ROC2_adcPulseTimeRaw":
                T_coin_pFADC_TREF_ROC2_adcPulseTimeRaw = e_tree.array("T.coin.pFADC_TREF_ROC2_adcPulseTimeRaw")
                treeDict.update({"T_coin_pFADC_TREF_ROC2_adcPulseTimeRaw" : T_coin_pFADC_TREF_ROC2_adcPulseTimeRaw})
            if branch == "T_coin_hFADC_TREF_ROC1_adcPulseTimeRaw":
                T_coin_hFADC_TREF_ROC1_adcPulseTimeRaw = e_tree.array("T.coin.hFADC_TREF_ROC1_adcPulseTimeRaw")
                treeDict.update({"T_coin_hFADC_TREF_ROC1_adcPulseTimeRaw" : T_coin_hFADC_TREF_ROC1_adcPulseTimeRaw})
            if branch == "T_coin_pEDTM_tdcTimeRaw":
                T_coin_pEDTM_tdcTimeRaw = e_tree.array("T.coin.pEDTM_tdcTimeRaw")
                treeDict.update({"T_coin_pEDTM_tdcTimeRaw" : T_coin_pEDTM_tdcTimeRaw})
            if branch == "T_coin_pEDTM_tdcTime":
                T_coin_pEDTM_tdcTime = e_tree.array("T.coin.pEDTM_tdcTime")
                treeDict.update({"T_coin_pEDTM_tdcTime" : T_coin_pEDTM_tdcTime})
                
            # Misc quantities
            if branch == "RFFreq":
                RFFreq = e_tree.array("MOFC1FREQ")  
                treeDict.update({"RFFreq" : RFFreq})
            if branch == "RFFreqDiff":
                RFFreqDiff = e_tree.array("MOFC1DELTA")
                treeDict.update({"RFFreqDiff" : RFFreqDiff})
            if branch == "EvtType":
                EvtType = e_tree.array("fEvtHdr.fEvtType")
                treeDict.update({"EvtType" : EvtType})
                
        #################################################################################################################
            
        # For better explaination of the methods below use the Help class defined above
        cutNames = []        
        cutVals = []
        if self.cuts != None:
            # read in cuts file and makes dictionary
            importDict = SetCuts(self.CURRENT_ENV).importDict(self.cuts,self.cut_f,self.runNum,self.DEBUG)
            for i,cut in enumerate(self.cuts):
                # Converts the dictionary to a list of strings that need to be evaluated and converted
                # into a boolean list
                x = SetCuts(self.CURRENT_ENV,importDict).booleanDict(cut)
                print("\n%s" % cut)
                print(x, "\n")
                # Saves string names of cuts and their values
                cutNames.append(cut)
                cutVals.append(x)
                if i == 0:
                    inputDict = {}
                # Redefines the dictionary to be reimplemented below
                cutDict = SetCuts(self.CURRENT_ENV,importDict).readDict(cut,inputDict)
                for j,val in enumerate(x):
                    try:
                        # Evaluates the list of strings which converts them to a list of boolean values
                        # corresponding to the cuts applied
                        cutDict = SetCuts(self.CURRENT_ENV,importDict).evalDict(cut,eval(x[j]),cutDict)
                        # This is the cython defined version, slightly faster but 
                        # requires it to be compiled fist
                        #cutDict = evalDict(cut,eval(x[j]),cutDict)
                    except NameError:
                        if "pid" in x[j]:
                            err_dir = self.UTILPATH+"/DB/PARAM/PID_Parameters.csv"
                        if "track" in x[j]:
                            err_dir = self.UTILPATH+"/DB/PARAM/Tracking_Parameters.csv"
                        if "accept" in x[j]:
                            err_dir = self.UTILPATH+"/DB/PARAM/Acceptance_Parameters.csv"
                        if "coin_time" in x[j]:
                            err_dir = self.UTILPATH+"/DB/PARAM/Timing_Parameters.csv"
                        if "current" in x[j]:
                            err_dir = self.UTILPATH+"/DB/PARAM/Current_Parameters.csv"
                        if "misc" in x[j]:
                            err_dir = self.UTILPATH+"/DB/PARAM/Misc_Parameters.csv"
                        raise InvalidEntry('''
                        ======================================================================
                          ERROR: %s invalid.

                          Improperly defined cut at... 
                          %s
                        ----------------------------------------------------------------------
                          Check that run number %s is properly defined in...
                          %s
                        ======================================================================
                        ''' % (cut,x[j],self.runNum,err_dir))
            strDict = dict(zip(cutNames,cutVals))

            return [SetCuts(self.CURRENT_ENV,cutDict),treeDict,strDict]
        else:
            return [SetCuts(self.CURRENT_ENV),treeDict,None]
        
    def csv2root(inputDict,rootName):
        '''
        csv2root(inputDict,rootName)
                 |         |
                 |         --> rootName: Output root file name
                 --> inputDict: Input dictionary with csv data to be converted to root

        ----------------------------------------------------------------------------------------------
        Converts csv file to root file. Save arrays,lists,etc. from csv to root file as histograms
        '''
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

class Equations():
    '''        
    Equations()

    ----------------------------------------------------------------------------------------------
    
    This class stores a variety of equations often used in the KaonLT analysis procedure
    '''

    def missmass():
        '''
        missmass()

        ----------------------------------------------------------------------------------------------

        Define missing mass calculation. !!! Not currently implimented !!!
        '''
        print("missmass")

class Misc():
    '''
    Misc()

    ----------------------------------------------------------------------------------------------

    Current functions...
            - progressBar

    ----------------------------------------------------------------------------------------------

    Class of miscellaneous methods
    '''
    
    def progressBar(value, endvalue, bar_length=50):
        '''
        progressBar(value, endvalue, bar_length=50)
                    |      |         |
                    |      |         --> bar_length: Length of bar to output to terminal (default = 50)
                    |      --> endvalue: End of loop value - 1
                    --> value: Iteration value
                        
        ----------------------------------------------------------------------------------------------

        A simple progress bar to use in loops
        '''

        percent = float(value) / endvalue
        arrow = '=' * int(round(percent * bar_length)-1) + '>'
        spaces = ' ' * (bar_length - len(arrow))
        if percent == 1:
            endl = '\n'
        else:
            endl = ''

        sys.stdout.write(" \r[{0}] {1}%\r{2}".format(arrow + spaces, round(percent * 100), endl))
        sys.stdout.flush()

    @contextmanager
    def suppress_stdout():
        '''
        suppress_stdout()

        ----------------------------------------------------------------------------------------------

        Suppresses python output. Use in a with statement and everything within will be suppressed
        '''
        with open(os.devnull, "w") as devnull:
            old_stdout = sys.stdout
            sys.stdout = devnull
            try:  
                yield
            finally:
                sys.stdout = old_stdout

    def test_cpp():
        print('')
