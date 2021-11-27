#! /usr/bin/python

#
# Description: Just calls help functions for various methods to help users
# ================================================================
# Time-stamp: "2021-11-18 05:18:48 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import os,glob

from .ltsep import Root, Equations, Misc
from .cut import SetCuts
from .pathing import SetPath

class Help():
    '''
    Help()

    ----------------------------------------------------------------------------------------------

    Class that is used to help users setup their scripts and get information on various functions 
    used throughout the ltsep package
    '''
        
    def info(func):
        '''
        info(func)
             |
             --> func: Input class/function to call help()

        ----------------------------------------------------------------------------------------------

        Calls help() to get description of function
        '''
        help(func)
    
    def getDoc(func):
        '''
        getDoc(func)

        ----------------------------------------------------------------------------------------------

        Decorator that allows docstring to be used inside the function being called
        '''
        def wrapper(*args, **kwargs):
            return func(func, *args, **kwargs)
        return wrapper

    @getDoc
    def path_setup(path_setup):
        '''
        ----------------------------------------------------------------------------------------------
        For the pathing you do not need to define all of the keys in the dictionary (like below),
        rather choose which paths are being used in your specific string. Make sure all references
        to UTIL_PION or UTIL_KAONLT are defined using UTILPATH (or any other useful path listed below)
        ----------------------------------------------------------------------------------------------

        ################################################################################################################################################
        \'''
        ltsep package import and pathing definitions
        \'''

        import os
        import ltsep as lt

        # os.path.realpath(__file__) is your current directory path
        # This will grab the pathing for these variables based off the files in PATH_TO_DIR
        HCANAPATH=lt.SetPath(os.path.realpath(__file__)).getPath("HCANAPATH")
        REPLAYPATH=lt.SetPath(os.path.realpath(__file__)).getPath("REPLAYPATH")
        UTILPATH=lt.SetPath(os.path.realpath(__file__)).getPath("UTILPATH")
        PACKAGEPATH=lt.SetPath(os.path.realpath(__file__)).getPath("PACKAGEPATH")
        OUTPATH=lt.SetPath(os.path.realpath(__file__)).getPath("OUTPATH")
        ROOTPATH=lt.SetPath(os.path.realpath(__file__)).getPath("ROOTPATH")
        REPORTPATH=lt.SetPath(os.path.realpath(__file__)).getPath("REPORTPATH")
        CUTPATH=lt.SetPath(os.path.realpath(__file__)).getPath("CUTPATH")
        PARAMPATH=lt.SetPath(os.path.realpath(__file__)).getPath("PARAMPATH")
        SCRIPTPATH=lt.SetPath(os.path.realpath(__file__)).getPath("SCRIPTPATH")
        ANATYPE=lt.SetPath(os.path.realpath(__file__)).getPath("ANATYPE")
        USER=lt.SetPath(os.path.realpath(__file__)).getPath("USER")
        HOST=lt.SetPath(os.path.realpath(__file__)).getPath("HOST")

        ################################################################################################################################################

        print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER, HOST, REPLAYPATH))

        ################################################################################################################################################
        \'''
        Check that root/output paths and files exist for use
        \'''

        # Construct the name of the rootfile based upon the info we provided
        OUTPATH = UTILPATH+"/OUTPUT/Analysis/%sLT" % ANATYPE        # Output folder location
        rootName = UTILPATH+"/ROOTfiles/Analysis/Lumi/%s_%s_%s.root" % (ROOTPrefix,runNum,MaxEvent)     # Input file location and variables taking
        print ("Attempting to process %s" %(rootName))
        lt.SetPath(os.path.realpath(__file__)).checkDir(OUTPATH)
        lt.SetPath(os.path.realpath(__file__)).checkFile(rootName)
        print("Output path checks out, outputting to %s" % (OUTPATH))
        '''
        print(path_setup.__doc__)

    @getDoc
    def cut_setup(cut_setup):
        '''
        ----------------------------------------------------------------------------------------------
        Make sure you have the following in your script...
        ----------------------------------------------------------------------------------------------

        import uproot as up
        import ltsep as lt

        # ----> Add pathing variables as well (see lt.Help.path_setup() for more info)

        # ---> Add r = klt.Root() for converting back to root files (see lt.Help.info(lt.Root) for more info)

        # Convert root leaf to array with uproot
        # Array name must match what is defined in DB/CUTS/general/
        leaf_name  = tree.array("leaf.name") # The periods are replaced with underscores

        ################################################################################################################################################
        \'''
        Define and set up cuts
        \'''

        fout = "<path_to_run_type_cut>"

        cuts = ["runTypeCut1","runTypeCut2",<etc>,...]

        def make_cutDict(cuts,fout,runNum,CURRENT_ENV,DEBUG=False):
            \'''
            This method calls several methods in kaonlt package. It is required to create properly formated
            dictionaries. The evaluation must be in the analysis script because the analysis variables (i.e. the
            leaves of interest) are not defined in the kaonlt package. This makes the system more flexible
            overall, but a bit more cumbersome in the analysis script. Perhaps one day a better solution will be
            implimented.
            \'''

            # read in cuts file and make dictionary
            importDict = lt.SetCuts(CURRENT_ENV).importDict(cuts,fout,runNum,DEBUG=DEBUG)
            for i,cut in enumerate(cuts):
                x = lt.SetCuts(CURRENT_ENV,importDict).booleanDict(cut)
                print("\\n%s" % cut)
                print(x, "\\n")
                if i == 0:
                    inputDict = {}
                cutDict = lt.SetCuts(CURRENT_ENV,importDict).readDict(cut,inputDict)
                for j,val in enumerate(x):
                    cutDict = lt.SetCuts(CURRENT_ENV,importDict).evalDict(cut,eval(x[j]),cutDict)
            return lt.SetCuts(CURRENT_ENV,cutDict)

        c = make_cutDict(cuts,fout,runNum,os.path.realpath(__file__))

        # ---> If multple run type files are required then define a new run type file altogether. Do not try to 
        # chain run type files. It can be done, but is computationally wasteful and pointless.

        # To apply cuts to array...
        c.add_cut(array,"runTypeCut")
        '''
        print(cut_setup.__doc__)

    def searchPathFile(CURRENT_ENV):
        '''
        searchPathFile(CURRENT_ENV)
                       |
                       --> CURRENT_ENV: Input current enviroment path

        ----------------------------------------------------------------------------------------------
        Outputs the current enviroment file and its contents which establish script pathing
        '''

        CURRENT_ENV = CURRENT_ENV.replace(os.getlogin(),"${USER}") # Replace username with general ${USER} so it can be used broadly
        CURRENT_ENV = CURRENT_ENV.split("/UTIL_",1)[0] # Redefine path to hallc_replay_lt (if in replay env)
        CURRENT_ENV = CURRENT_ENV.split("/cut.py",1)[0] # Redefine path to ltsep (if in package env)

        # Grab location of package (either in .local or in the UTIL dir)
        PACKAGE_ENV = os.path.dirname(os.path.realpath(__file__))

        # Grab username and hostname
        USER = os.getlogin()
        HOST = os.uname()[1]

        # Setup path to pathing files (see PATH_TO_DIR)
        path_check = "{}/PATH_TO_DIR".format(PACKAGE_ENV)
        
        # Check through all pathing files for the one that matches the current working enviroment
        for fname in glob.glob(path_check+"/*.path"):
            with open(fname) as f:
                search = f.read()
            if USER == "cdaq":
                if PACKAGE_ENV in search:
                    PATHFILE = fname
            else:
                if CURRENT_ENV in search:
                    PATHFILE = fname

        print("\t----------------------------------------------------------------------------------------------")
        print("\tThe current enviroment path file used is...\n\t{}".format(PATHFILE))
        print("\t----------------------------------------------------------------------------------------------\n")
        with open(PATHFILE) as f:
            for line in f:
                print("\t",line)
        print("\n\n")
