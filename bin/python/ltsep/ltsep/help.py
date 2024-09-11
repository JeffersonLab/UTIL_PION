#! /usr/bin/python

#
# Description: Just calls help functions for various methods to help users
# ================================================================
# Time-stamp: "2022-06-30 02:35:57 trottar"
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
        ##############################################################################################
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
        ANATYPE=lt.ANATYPE
        USER=lt.USER
        HOST=lt.HOST
        # Note the OUTPATH is not defined unless RunType argument is given, see below

        # If you wish to explicitly define root branches then do the following...
        import uproot as up
        tree = up.open("<ROOT_FILE_NAME>")["<ROOT_TREE_NAME>"]
        # Convert root leaf to array with uproot
        branch_name  = tree.array("<ROOT_BRANCH_NAME>") # The periods are replaced with underscores

        ##############################################################################################
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
        ANATYPE=lt.ANATYPE
        USER=lt.USER
        HOST=lt.HOST
        OUTPATH=lt.OUTPATH

        # If you wish to explicitly define root branches then do the following...
        import uproot as up
        tree = up.open("<ROOT_FILE_NAME>")["<ROOT_TREE_NAME>"]
        # Convert root leaf to array with uproot
        branch_name  = tree.array("<ROOT_BRANCH_NAME>") # The periods are replaced with underscores

        ##############################################################################################
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

        ##############################################################################################
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

        ##############################################################################################
        \'''
        Define bash dynamic pathing
        \'''

        # Runs script in the ltsep python package that grabs current path enviroment
        if [[ ${HOSTNAME} = *"cdaq"* ]]; then
            PATHFILE_INFO=`python3 /home/cdaq/pionLT-2021/hallc_replay_lt/UTIL_PION/bin/python/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
        elif [[ "${HOSTNAME}" = *"farm"* ]]; then
            PATHFILE_INFO=`python3 /u/home/${USER}/.local/lib/python3.4/site-packages/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
        fi

        # Split the string we get to individual variables, easier for printing and use later
        VOLATILEPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f1` # Cut the string on , delimitter, select field (f) 1, set variable to output of command
        ANALYSISPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f2`
        HCANAPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f3`
        REPLAYPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f4`
        UTILPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f5`
        PACKAGEPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f6`
        OUTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f7`
        ROOTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f8`
        REPORTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f9`
        CUTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f10`
        PARAMPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f11`
        SCRIPTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f12`
        ANATYPE=`echo ${PATHFILE_INFO} | cut -d ','  -f13`
        USER=`echo ${PATHFILE_INFO} | cut -d ','  -f14`
        HOST=`echo ${PATHFILE_INFO} | cut -d ','  -f15`
        SIMCPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f16`

        ----------------------------------------------------------------------------------------------

        '''
        print(path_setup.__doc__)

    @getDoc
    def cut_setup(cut_setup):
        '''
        ----------------------------------------------------------------------------------------------
        Make sure you have the following in your script...
        ----------------------------------------------------------------------------------------------
        ##############################################################################################
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
        ##############################################################################################
        ----------------------------------------------------------------------------------------------

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
