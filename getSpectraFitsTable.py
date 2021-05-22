
#########################################
###      create `specResultsTable`    ###
#########################################

"""
22nd May 2021:
---
Have completed much of the pending work regarding Spectral Analysis
of ASASSN-14li. This Python routine is to read the `specModelParams.json`
files for each of the obsIDs and create a `specResultsTable` using them.
"""

import os
import json
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from astropy.table import Table

from epicObj import epicObj
from plotAnal import plotAnal


##-- the main function --##
def main (args):
    
    #-- create an object of class spectra --#
    obj = epicObj(ra=args.ra, dec=args.dec, workdir=args.workdir)
    
    #-- get obsIDs and the objName --#
    if args.objName == None:
        try:
            obj.findObsIDs()
        except KeyError:
            print('\nNo obsIDs found for the given location. Confirm that you are connected to the Internet!')
    else:
        obj.objName = args.objName

    #-- check if obsID has been passed --#
    if args.obsIDs!=None:
        obj.obsIDs = args.obsIDs
        print('Using the obsIDs passed.')

    return print('\nSuccessfully created the specResultsTable!')


if __name__=="__main__":

    description = 'Python routine to create specResultsTable using specModelParams.json'
    parser = argparse.ArgumentParser(description=description)

    #-- general arguments --#
    parser.add_argument('--ra', action='store', type=float, default=192.0625, \
                        help='right ascension of the object. (default:%(default)s, ASASSN-14li)')
    parser.add_argument('--dec', action='store', type=float, default=17.7739, \
                        help='declination of the object. (default:%(default)s, ASASSN-14li)')    
    parser.add_argument('--workdir', action='store', type=str, default='/media/suyog/DATA/xmm_obs', \
                        help='directory where obsid folders will be stored. (default:%(default)s)')
    parser.add_argument('--obsIDs', nargs='+', action='store', default=None, #['0831790201'], \
                        help='''obsIDs for which some specific function has to executed. 
                                Valid only when --method argument is specified. (default:%(default)s)''')
    parser.add_argument('--objName', action='store', default=None, \
                        help='name of the obj used to locate the xmmObj pickle file.')

    args = parser.parse_args()

    main(args)







#################### End of Program #########################
#############################################################
