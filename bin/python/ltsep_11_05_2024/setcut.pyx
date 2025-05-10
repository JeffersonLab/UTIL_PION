#! /usr/bin/python
#
# Description:
# ================================================================
# Time-stamp: "2022-06-15 12:57:08 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

#print('Hello World!')

def booleanDict(inputDict,cut):
    '''
    booleanDict(self,cut)
                     |
                     --> cut: Input cut to convert to boolean

    ----------------------------------------------------------------------------------------------

    Create a boolean dictionary for cut by converting string to array of pass/no pass cut.
    '''

    subDict = inputDict[cut]
    subDict = subDict.split(",")
    cut_arr = [evt for evt in subDict]
    cut_arr = list(filter(None,cut_arr)) # Filter out blank string values (ie where cuts were removed)
    return cut_arr

def evalDict(cut,eval_xi,inputDict):
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

def apply_cut(inputDict,arr, cut):
    '''
    apply_cut(self,arr, cut)
                 |    |
                 |    --> cut: Run type cut name to impliment to array 
                 --> arr: Input array to be cut

    ----------------------------------------------------------------------------------------------

    Creates the string of cuts to later be evaluated and then converted to boolean array.
    '''

    applycut = "arr["
    subDict = inputDict[cut]
    for i,(key,val) in enumerate(subDict.items()):
        if i == len(subDict)-1:
            applycut += 'self.cut("%s","%s")]' % (key,cut)
        else:
            applycut += 'self.cut("%s","%s") & ' % (key,cut)
    return applycut
