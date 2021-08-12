'''
Sttering script to run SIMC automatically and generate 
an output ROOTfile for further analysis.

Author : C. Yero
Date   : May 17, 2021
email  : cyero@jlab.org
'''

#import numpy as np    
import os             # contains the system function that can execute shell commands 

# input file name containng list of kinematic files for simulation
input_list_name = "./input_list.txt"
# read input kin. file list
input_list = open(input_list_name, 'r')

#loop over every kin. file in the list
for ikin in input_list:

    # ignore commented lines
    if(ikin[0]=='#'):
        continue
    else:
        #remove blanks
        ikin = ikin.strip() 

        #define shell script command
        run_simc_cmd = "./run_simc %s " % (ikin)
        print(run_simc_cmd)

        #execute shell script command to run simulation for ikin
        os.system(run_simc_cmd)

        #define command to convert ASCII file to ROOTfile
        run_fmake_cmd = "./run_simc_tree %s " % (ikin) 
        print(run_fmake_cmd)

        #execute command to convert to ROOTfile
        os.system(run_fmake_cmd)
