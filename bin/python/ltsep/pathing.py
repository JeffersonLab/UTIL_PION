#! /usr/bin/python
# Description: This package will perform many tasks required for l-t separation physics analysis 
# Analysis script required dynamically defining pathing.
# ================================================================
# Time-stamp: "2021-11-05 07:50:04 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
from pathlib import Path
import sys,os,glob

class SetPath():
    '''
    SetPath()

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

    ----------------------------------------------------------------------------------------------

    Class that sets the pathing for scripts as well as finds if dir, symlink, or file exists
    '''

    def __init__(self, CURRENT_ENV):
        '''
        __init__(self,CURRENT_ENV):
                      |
                      --> CURRENT_ENV: Input current enviroment path

        ----------------------------------------------------------------------------------------------

        Constructor of class takes the current enviroment path as input
        '''

        CURRENT_ENV = CURRENT_ENV.split("/UTIL_",1)[0] # Redefine path to hallc_replay_lt (if in replay env)
        CURRENT_ENV = CURRENT_ENV.split("/cut.py",1)[0] # Redefine path to ltsep (if in package env)
        self.CURRENT_ENV = CURRENT_ENV

    def __str__(self):
        '''
        __str__(self)

        ----------------------------------------------------------------------------------------------

        String representation of class if called as string (eg print(SetPath))
        '''

        return "CURRENT_ENV : {self.CURRENT_ENV}"

    def __repr__(self):
        '''
        __repr__(self)

        ----------------------------------------------------------------------------------------------

        String representation of class if called as is (eg SetPath)
        '''

        return "SetCuts({self.CURRENT_ENV})"

    def getPath(self,inp_dir,DEBUG=False):
        '''
        getPath(self,inp_dir,DEBUG=False)
                     |       |
                     |       --> DEBUG: Debugging flag
                     ----------> inp_dir: Key to dictionary

        ----------------------------------------------------------------------------------------------

        Get path of working directory and set up dictionary with pathing strings
        '''

        # Grab location of package (either in .local or in the UTIL dir)
        PACKAGE_ENV = os.path.dirname(os.path.realpath(__file__))

        # If python package in user dir, then get username from string
        # This is required for batch jobs to run properly (in practice amounts to the same as os.getlogin())
        if "local" in PACKAGE_ENV:
            USER = PACKAGE_ENV.split("/.local")[0]
            USER = USER.split("home/")[1]
        else:
            USER = os.getlogin()

        # Grab username and hostname
        HOST = os.uname()[1]

        # Replace username with general ${USER} so it can be used broadly
        self.CURRENT_ENV = self.CURRENT_ENV.replace(USER,"${USER}") 

        if DEBUG==True:
            print("CURRENT_ENV ",self.CURRENT_ENV)

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
                if self.CURRENT_ENV in search:
                        PATHFILE = fname
        try:
            PATHFILE
        except NameError:
            print("ERROR: PATHFILE not defined. Invalid enviroment...\n\t{}".format(self.CURRENT_ENV))
            sys.exit(1)

        # Open pathing file then create a pathing dictionary based off the contents
        inp_path = open(PATHFILE)
        pathDict = {}
        for line in inp_path:
            line  = line.split("=")
            # Create dictionary of pathing. Assuring that ${USER} is replaced with proper user names
            pathDict.update({line[0].strip().strip("\n") : line[1].strip().strip("\n").replace("${USER}",USER)})
        # Include username and hostname to dictionary
        pathDict.update({"USER" : USER})
        pathDict.update({"HOST" : HOST})
        if DEBUG==True:
            print("pathDict ",pathDict)
        inp_path.close()

        return pathDict[inp_dir]

    def checkDir(self,inp_dir):
        '''
        checkDir(self,inp_dir)
                      |
                      --> inp_dir: Input dir/symlink to check

        ----------------------------------------------------------------------------------------------

        Check if directory and/or symlink exist
        '''

        if os.path.exists(inp_dir):
            if os.path.islink(inp_dir):
                pass
            elif os.path.isdir(inp_dir):
                pass
            else:
                print ("{} exists but is not a directory or sym link, check your directory/link and try again".format(inp_dir))
                sys.exit(2)
        else:
            print("Output path not found, please make a sym link or directory called OUTPUT in {} to store output").format(UTILPATH.replace(REPLAYPATH+"/",""))
            sys.exit(3)

    def checkFile(self,inp_file):
        '''
        checkFile(self,inp_file)        
                       |
                       --> inp_file: Input file to check

        ----------------------------------------------------------------------------------------------

        Check if file exists
        '''

        if os.path.isfile(inp_file):
            print ("{} exists, processing".format(inp_file))
        else:
            print ("{} not found - do you have the correct sym link/folder set up?".format(inp_file))
            sys.exit(4)
