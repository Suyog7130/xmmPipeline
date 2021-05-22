
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

from datetime import datetime
from astropy.table import Table

from epicObj import epicObj
from plotAnal import plotAnal


##-- the main function --##
def main (args):
    
    #-- create an object of class spectra --#
    obj = epicObj(ra=args.ra, dec=args.dec, workdir=args.workdir, \
                  sas_dir=None, headas=None, sas_ccfpath=None)

    #-- check if obsID has been passed --#
    if args.obsIDs != None:
        obsIDs = args.obsIDs
        print('Using the obsIDs passed.')
    else:
        obj.findObsIDs()
        obsIDs = obj.obsIDs

    #-- get sorted obsIDs --#
    df = obj.sortObsIDs()
    df.time = [t.date() for t in df.time]

    #-- the models to use --#
    models = ['tbabs*clumin*zashift*(bbobyrad+bbodyrad)', 'tbabs*clumin*zashift*bbodyrad', \
              'tbabs*clumin*zashift*(bbodyrad+powerlaw)', 'tbabs*clumin*zashift*powerlaw', \
              'tbabs*cflux*zashift*(bbodyrad+bbodyrad)', 'tbabs*cflux*zashift*bbodyrad', \
              'tbabs*cflux*zashift*(bbodyrad+powerlaw)', 'tbabs*cflux*zashift*powerlaw', \
              'tbabs*clumin*zashift*(bbodyrad+powerlaw+diskbb)']
    if not args.model == 'all':
        models = [args.model]
    
    #-- iterate for all the models --#
    #for model in models:
        

    print(df)
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

    #-- fitSpectra arguments --#
    parser.add_argument('--model', action='store', default='tbabs*zashift*(powerlaw)', \
                        help='name of the model to be used. (default:%(default)s)')

    args = parser.parse_args()

    main(args)







#################### End of Program #########################
#############################################################
