
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

23rd May 2021:
---
The primary ``for`` loop should be for model since each model will have a
separate ``specResultsTable``.

18th June 2021:
---
Finishing up the manual specModel fitting.
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

pd.set_option('display.expand_frame_repr', False)


##-- the main function --##
def main (args):
    verbose = args.verbose
    
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
    df.time = [str(t.date()) for t in df.time]
    df.columns = ['obsIDs', 'date']
    sortedObsIDs = df.obsIDs.tolist()
    df = df.set_index('obsIDs')
    #df['date'].astype(str)

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
        if verbose:
            print(f'\nCreating specResultsTable for model\n\t{model}')
    
        #-- iterate for each obsIDs --#
        for obsID in sortedObsIDs:
            if verbose:
                print(f'\nLoading specModelParams.json for {obsID}')
            workdir = args.workdir +'/'+ obsID +'/work'
            
            #-- load the `specModelParams.json` file --#
            fname = workdir+'/specModelParams.json'
            if os.path.isfile(fname):
                resultDict = json.load( open(fname, 'r') )
                #print(resultDict)
            else:
                if verbose:
                    printErrorMessage('file not found!')
                    print(f'\nRow for {obsID} will not be added.')
                continue

            #-- check for presence of the model --#
            modelDict = resultDict.get(model, None)
            if modelDict is None:
                if verbose:
                    print('\nPresent model is not saved in the JSON file.')
                continue

            #-- add columns for useful model parameters --#
            df = df.fillna('None')
            tableCols = ['kT', 'norm', 'PhoIndex', 'norm', \
                         'lg10Lum', 'lg10Flux', 'pnLumin', 'pnFlux', 'chiSq', 'dof']
            mComps = [k for k in modelDict.keys() if k != 'results']

            for comp in mComps:
                mParams = modelDict[comp].keys()
                for param in mParams:
                    if param in tableCols:
                        pCol = param
                        if param in df.columns.tolist() and df.loc[obsID, param] != 'None':
                            pCol = param + '_1'
                        pErrorCol = pCol + '_sigma'
                        df.loc[obsID, pCol] = modelDict[comp][param]['value']
                        df.loc[obsID, pErrorCol] = modelDict[comp][param]['sigma']

            #-- add result metric columns --#
            metricDict = modelDict['results']
            for metric in metricDict.keys():
                if metric in tableCols:
                    val = metricDict[metric]
                    if metric in ['pnLumin', 'mos1lumin', 'mos2lumin']:
                        val = val['eUnit_Ine44']
                    if metric in ['pnFlux', 'mos1Flux', 'mos2Flux']:
                        val = val['eUnit']
                    df.loc[obsID, metric] = val

        #-- if only some obsIDs have to be included --#
        if args.obsIDs is not None:
            toDrop = list( set(df.index.tolist()) - set(args.obsIDs) )
            df = df.drop(toDrop)
            if verbose:
                print('Using the obsIDs passed.')

        df = df.replace('None', 0.0)  #--convert empty values to float64.

        #-- move ``chiSq`` and ``dof`` columns to the end --#
        kTcols = ['kT', 'kT_sigma', 'norm', 'norm_sigma']
        phoCols = ['PhoIndex', 'PhoIndex_sigma', 'norm_1', 'norm_1_sigma']
        pho0cols = ['PhoIndex', 'PhoIndex_sigma', 'norm', 'norm_sigma']
        kT1cols = ['kT_1', 'kT_1_sigma', 'norm_1', 'norm_1_sigma']
        pho1cols = ['PhoIndex_1', 'PhoIndex_1_sigma', 'norm_1', 'norm_1_sigma']
        otherCols = ['lg10Lum', 'lg10Lum_sigma', 'lg10Flux', 'lg10Flux_sigma', \
                     'pnLumin', 'pnFlux', 'chiSq', 'dof']

        dfColsList = df.columns.tolist()
        if 'kT' in dfColsList and 'kT_1' in dfColsList:
            moveColToEndOrder = kTcols + kT1cols + otherCols
        elif 'kT' in dfColsList and 'PhoIndex' in dfColsList:
            moveColToEndOrder = kTcols + phoCols + otherCols
        else:
            moveColToEndOrder = otherCols

        if args.orderTableCols:
            for col in moveColToEndOrder:
                if col in dfColsList:
                    df.insert(len(df.columns)-1, col, df.pop(col))

        #-- add the ``redChiSq`` value column --#
        redChiSq = []
        for chiSq, dof in zip(df.chiSq.tolist(), df.dof.tolist()):
            if dof != 0:
                redChiSq.append(chiSq/dof)
            else:
                redChiSq.append(0)
        df['redChiSq'] = redChiSq

        print(f'\nThe final specResultsTable for model\n\t{model}\n')  
        print(df)   #, df.columns.tolist())
        print('\nSuccessfully created the specResultsTable!')

        df.to_csv('specResultsTable_'+model+'.csv')
        if verbose:
            print('Saved the specResultsTable to CSV.')

    return True


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
                        help='''obsIDs for which the \'specResultsTable\' has to obtained. 
                                (default:%(default)s)''')
    parser.add_argument('--objName', action='store', default=None, \
                        help='name of the obj used to locate the xmmObj pickle file.')

    parser.add_argument('-v', '--verbose', action='store_true', default=False, \
                        help='more messages shown on the terminal. (default:%(default)s)')
    parser.add_argument('--orderTableCols', action='store_true', default=False, \
                        help='order the specResultsTable columns. (default:%(default)s)')

    #-- fitSpectra arguments --#
    defaultModel = 'tbabs*clumin*zashift*(bbodyrad+bbodyrad)_1frozen'
    parser.add_argument('--model', action='store', default=defaultModel, \
                        help='name of the model to be used. (default:%(default)s)')

    args = parser.parse_args()

    main(args)







#################### End of Program #########################
#############################################################
