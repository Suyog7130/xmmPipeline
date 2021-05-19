
######################################
###      Reduce the XMM Data       ###
######################################

"""
10th May 2021:
---
Am making several changes to the code procedure arrangements. 
See the GitHub repo, the Notes on Google Docs and the documentation for more information.
For previous docstring comments, see earlier code files, namely `xmmPipeline.py`
"""

import os
import subprocess
import requests
import wget

import glob
import pickle
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from astropy.table import Table

from epicObj import epicObj
from epicObj import strToBool, printErrorMessage
from plotAnal import plotAnal


##-- the reduceData main function --##
def reduceData_method (args):

    #-- create an object of class xmmObj --#
    obj = epicObj(ra=args.ra, dec=args.dec, workdir=args.workdir, \
                  sas_dir=args.sas_dir, headas=args.headas, sas_ccfpath=args.sas_ccfpath)
                 
    #-- run the reduceData functions --#
    try:
        obj.findObsIDs()
    except KeyError:
        print('\nNo obsIDs found for the given location. Confirm that you are connected to the Internet!')
    obj.downloadData()
    obj.reduceEPICdata()
    
    return print('\nHurray! reduceData Method Successfully ran.') 
    
    
##-- the extractProds main function --##
def combineAndFind_method (args):

    if args.bkgCircRadius is not None and len(args.bkgCircRadius)==1:
        r = float(args.bkgCircRadius)
        bkgCircRadius = [r, r]
    else:
        bkgCircRadius = args.bkgCircRadius

    #-- create an object of class xmmObj --#
    obj = epicObj(ra=args.ra, dec=args.dec, workdir=args.workdir, \
                  sas_dir=args.sas_dir, headas=args.headas, sas_ccfpath=args.sas_ccfpath, \
                  srcCircRadius=args.srcCircRadius, bkgCircRadius=bkgCircRadius, dSrcThreshold=args.dSrcThreshold, \
                  bkgCircGap=args.bkgCircGap, esp_nsplinenodes=args.esp_nsplinenodes, \
                  gti_indiThreshold=args.gti_indiThreshold, gti_combThreshold=args.gti_combThreshold, \
                  lcBinSize=args.lcBinSize, binBkglc=args.binBkglc, \
                  saveFig=args.saveFig, showFig=args.showFig)
    
    #-- check if obsID has been given --#
    if args.obsIDs!=None:
        obj.obsIDs = args.obsIDs
    else:
        try:
            obj.findObsIDs()
        except KeyError:
            print('\nNo obsIDs found for the given location. Confirm that you are connected to the Internet!')
            
    #-- run the extract products functions --#
    obj.removeFlareBackground()
    obj.findSourceCCD()
    obj.extract_InstrumentalGTIs()
    obj.combineGTIs()
    obj.combineEPICdata()
    obj.find_otherSources()
    obj.getBackgroundCircles()
    obj.writeCCDcoordsPickle()

    print('\nHurray! extractProds Method Succesfully ran.')
    if len(obj.badObs)!=0:
        print('These obsIDs were excluded from analysis: ', obj.badObs)

    return True


##-- the extractProds main function --##
def extractProds_method (args):

    #-- create an object of class xmmObj --#
    obj = epicObj(ra=args.ra, dec=args.dec, workdir=args.workdir, \
                  sas_dir=args.sas_dir, headas=args.headas, sas_ccfpath=args.sas_ccfpath, \
                  lcBinSize=args.lcBinSize, binBkglc=args.binBkglc, \
                  saveFig=args.saveFig, showFig=args.showFig, ignorePileup=args.ignorePileup)
    
    #-- check if obsID has been given --#
    if args.obsIDs!=None:
        obj.obsIDs = args.obsIDs
    else:
        try:
            obj.findObsIDs()
        except KeyError:
            print('\nNo obsIDs found for the given location. Confirm that you are connected to the Internet!')

    if args.saveResults:
        obj.save_results()
        return True
            
    #-- run the extract products functions --#
    obj.readCCDcoordsPickle()
    obj.extract_srcBkg_eventLists()
    #obj._extract_sMode_srcPNlc(obsID=obj.obsIDs[0])
    obj.obtain_lightCurves()
    obj.save_results()

    print('\nHurray! extractProds Method Succesfully ran.')
    if len(obj.badObs)!=0:
        print('These obsIDs were excluded from analysis: ', obj.badObs)

    print('\nNow saving this whole epicObj to a Pickle file.')
    outfile = open(args.workdir+'/'+'output_epicObj.pickle', 'wb')
    pickle.dump(obj, outfile)
    outfile.close()
    print('Saved to Pickle file!')

    return True


##-------------------------------------------------------------------------------------------##

if __name__=="__main__":
    
    description = 'Download, Reduce and extract products from the XMM data of the object located at the given coordinates.'
    
    parser = argparse.ArgumentParser(description=description)   #--create a ArgumentParser object.

    #-- general arguments --#
    parser.add_argument('--ra', action='store', type=float, default=192.0625, \
                        help='right ascension of the object. (default:%(default)s, ASASSN-14li)')
    parser.add_argument('--dec', action='store', type=float, default=17.7739, \
                        help='declination of the object. (default:%(default)s, ASASSN-14li)')   
    parser.add_argument('--workdir', action='store', type=str, default='/media/suyog/DATA/xmm_obs', \
                        help='directory where obsid folders will be stored. (default:%(default)s)')
            
    #-- arguments for running specific functions --#            
    parser.add_argument('--method', action='store', type=str, #default='reduceData', \
                        choices=['reduceData', 'combineAndFind', 'extractProds'], \
                        help='the set of specific functions to be executed. Default method is set to None.')
    parser.add_argument('--obsIDs', nargs='+', action='store', default=None, #['0831790201'], \
                        help='''obsIDs for which some specific function has to executed. 
                                Valid only when --method argument is specified. (default:%(default)s)''')
    
    #-- location paths arguments --#
    SAS_DIR = '/usr/local/xmmsas_20201028_0905'  
    HEADAS = '/usr/local/heasoft-6.28/x86_64-pc-linux-gnu-libc2.27'
    SAS_CCFPATH = '/ccf'
    parser.add_argument('--sas_dir', action='store', type=str, default=SAS_DIR, \
                        help='SAS_DIR environment variable. (default:%(default)s)')
    parser.add_argument('--headas', action='store', type=str, default=HEADAS, \
                        help='HEADAS environment variable. (default:%(default)s)')
    parser.add_argument('--sas_ccfpath', action='store', type=str, default=SAS_CCFPATH, \
                        help='SAS_CCFPATH environment variable. (default:%(default)s)')

    #-- optional arguments for function parameters --#
    group = parser.add_argument_group('function parameters')
    group.add_argument('--srcCircRadius', action='store', default=None, \
                        help='radius of the source circle, arcsec. \
                              Default value is half of dSrcThreshold or the perpendicular distance to the nearest boundary, \
                              whichever is minimum.')
    group.add_argument('--bkgCircRadius', action='store', nargs='+', default=None, \
                        help='radii of the background circles, arcsec. \
                              Default values are r1=0.5*threshold and r2=r1-10.')
    group.add_argument('--dSrcThreshold', action='store', default=70, \
                        help='threshold distance from the source, arcsec. (default:%(default)s)')
    group.add_argument('--bkgCircGap', action='store', default=2.5, \
                        help='gap to have around the circles, in pixels. (default:%(default)s)')

    #group.add_argument('--edetectmode', action='store', choices=['individual', 'chain'], default='chain', \
    #                    help='run edetect_chain or individual edetect tasks?')
    group.add_argument('--esp_nsplinenodes', action='store', default=14, \
                        help='for esplinemap within edetect_chain, controls no. of sources detected. (default:%(default)s)')
    group.add_argument('--gti_indiThreshold', action='store', default=500, \
                        help='include individuals GTIs greater than this value, seconds. (default:%(default)s)')
    group.add_argument('--gti_combThreshold', action='store', default=1000, \
                        help='remove obsIDs having total GTI below this value from further analysis, seconds. \
                              Default value is twice gti_indiThreshold.')
    group.add_argument('--lcBinSize', action='store', default=25, 
                        help='binning size for the light curve, seconds. (default:%(default)s). \
                              Passing (linear) makes the srcBinSize equal to 2 times the max Source count rate, \
                              while the bkgBinSize takes the default value of 25 sec.')
    group.add_argument('--binBkglc', action='store', default='no', \
                        help='whether to bin background light curve or not? (default:%(default)s)')

    parser.add_argument('--saveFig', action='store_true', default=False, \
                        help='save the matplotlib plots or not? (default:%(default)s)')
    parser.add_argument('--showFig', action='store_true', default=False, \
                        help='show the matplotlib plots or not? (default:%(default)s)')
    parser.add_argument('--saveResults', action='store_true', default=False, \
                        help='run only save_results function. (default:%(default)s)')
    parser.add_argument('--ignorePileup', action='store_true', defalt=False, \
                        help='ignore Pile-up in Piled-up obsIDs. (default:%(default)s)')
    
    #-- parse the arguments --#
    args = parser.parse_args()   #--parse all the arguments.
    #print(args)
    
    #-- select the pipeline method --#
    if 'func' in args:
        args.func(args)
    elif args.method == 'reduceData':
        reduceData_method(args)
    elif args.method == 'combineAndFind':
        combineAndFind_method(args)
    elif args.method == 'extractProds':
        extractProds_method(args)
    elif args.saveResults:
        extractProds_method(args)
    else:
        printErrorMessage('Please give which method to proceed with!')
        parser.print_help()

   



#################### End of Program #########################
#############################################################



