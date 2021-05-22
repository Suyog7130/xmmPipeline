
#########################################
###      create `specResultsTable`    ###
#########################################

"""
22nd May 2021:
---
Have completed much of the pending work regarding Spectral Analysis
of ASASSN-14li. This Python routine is to read the `specModelParams.json`
files for each of the obsIDs and create a `specResultsTable` using them.

I don't the obsIDs thing will be valid here.
Why would one not wanna include some obsIDs in the final table?
Anyway.
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
from epicObj import strToBool, printErrorMessage
from plotAnal import plotAnal


##-- the main function --##
def main (args):
    
    #-- create an object of class spectra --#
    obj = epicObj(ra=args.ra, dec=args.dec, workdir=args.workdir, \
                  sas_dir=None, headas=None, sas_ccfpath=None)

    #-- check if objName has been passed --#
    if args.objName is None:
        obj.findObsIDs()
        df = obj.sortObsIDs()
    else:
        df = obj.sortObsIDs(objName=args.objName)

    #-- get sorted obsIDs --#
    df.time = [t.date() for t in df.time]
    df.columns = ['obsIDs', 'date']
    sortedObsIDs = df.obsIDs.tolist()
    df = df.set_index('obsIDs')

    #-- the models to use --#
    models = ['tbabs*clumin*zashift*(bbobyrad+bbodyrad)', 'tbabs*clumin*zashift*bbodyrad', \
              'tbabs*clumin*zashift*(bbodyrad+powerlaw)', 'tbabs*clumin*zashift*powerlaw', \
              'tbabs*cflux*zashift*(bbodyrad+bbodyrad)', 'tbabs*cflux*zashift*bbodyrad', \
              'tbabs*cflux*zashift*(bbodyrad+powerlaw)', 'tbabs*cflux*zashift*powerlaw', \
              'tbabs*clumin*zashift*(bbodyrad+powerlaw+diskbb)']
    if not args.model == 'all':
        models = [args.model]
    
    #-- iterate for all the models --#
    for model in models:
        print(f'\nCreating specResultsTable for model\n{model}')

        #-- iterate for each obsIDs --#
        for obsID in sortedObsIDs:
            print(f'\nLoading specModelParams.json for {obsID}')
            workdir = args.workdir +'/'+ obsID +'/work'
            
            #-- load the `specModelParams.json` file --#
            fname = workdir+'/specModelParams.json'
            if os.path.isfile(fname):
                resultDict = json.load( open(fname, 'r') )
                #print(resultDict)
            else:
                printErrorMessage('file not found!')
                print(f'\nRow for {obsID} will not be added.')
                continue

            #-- check for presence of the model --#
            modelDict = resultDict.get(model, None)
            if modelDict is None:
                print('\nPresent model is not saved in the JSON file.')
                continue

            #-- add columns for useful model parameters --#
            mComps = [k for k in modelDict.keys() if k != 'results']
            for comp in mComps:
                mParams = modelDict[comp].keys()
                for param in mParams:
                    if param in ['kT', 'norm', 'lg10Lum']:
                        if param in df.columns.tolist():
                            pCol = param + '_1'
                        else:
                            pCol = param
                        df.loc[obsID, pCol] = modelDict[comp][param]['value']

            #-- add result metric columns --#
            

    #-- if only some obsIDs have to be included --#
    if args.obsIDs is not None:
        toDrop = list( set(df.index.tolist()) - set(args.obsIDs) )
        df = df.drop(toDrop)
        print('Using the obsIDs passed.')

    print(df)  #, df.columns.tolist())
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
    defaultModel = 'tbabs*clumin*zashift*(bbodyrad+bbodyrad)_1frozen'
    parser.add_argument('--model', action='store', default=defaultModel, \
                        help='name of the model to be used. (default:%(default)s)')

    args = parser.parse_args()

    main(args)







#################### End of Program #########################
#############################################################
