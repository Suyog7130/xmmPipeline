
####################################
###      EPIC Pile-up Issue      ###
####################################

"""
08 May 2021:
---
For clearing the EPIC Pile-up in early time obsIDs.
See: `https://www.cosmos.esa.int/web/xmm-newton/sas-thread-epatplot`

This would be used to correct for the Pile-up seperately for individual
obsIDs. Thus, the `workdir` is changed to `workdir`+'/work' folder.
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

from xmmPipeline import xmmObj
from spectra import spectra
from plotAnal import plotAnal


##-- function to convert yes/no to bool --##
def strToBool (s):
    if type(s)==bool:
        return s
    elif s in ['yes', 'y', 'true', 'True', 'Y', 'YES', 'TRUE']:
        return True
    elif s in ['no', 'n', 'false', 'False', 'N', 'NO', 'FALSE']:
        return False
    else:
        return print('\nPlease give bool values as yes/no.')

##-- print error message --##
def printErrorMessage (message):
    width = len(str(message))+4
    message = str(message).center(width, ' ')
    print('\n\t\t'+'*'*(width+4))
    print(f'\t\t**{message}**')
    print('\t\t'+'*'*(width+4))


class epicPileup (metaclass=spectra):

    def checkPileUp (self, workdir, obsID, saveFig=True, showFig=True):

        #-- set the environment variables --#
        os.environ['SAS_DIR'] = self.sas_dir
        os.environ['HEADAS'] = self.headas
        os.environ['SAS_CCFPATH'] = self.sas_ccfpath
        
        #-- iterating for all the obsIDs --#      
        for obsID in self.obsIDs:
            print('\nExtracting Image Mode Spectra for obsID {}.'.format(obsID))
            workdir = self.workdir+'/'+obsID+'/work'
            
            #-- get the location parameters --#
            locParams = self.sourceLoc[obsID]
            srcX, srcY, srcR = str(locParams['x']), str(locParams['y']), str(locParams['r'])
            
            bLocParams = self.backgroundLoc[obsID]
            Bx1, By1, Br1 = str(bLocParams['Bx1']), str(bLocParams['By1']), str(bLocParams['Br1'])
            Bx2, By2, Br2 = str(bLocParams['Bx2']), str(bLocParams['By2']), str(bLocParams['Br2'])
            #print(srcX, srcY, srcR, '\n', Bx1, By1, Br1, '\n', Bx2, By2, Br2)

            #-- get the CCD numbers --#
            pnCCD = str(self.sourceCCDs[obsID]['PN'])
            mos1CCD = str(self.sourceCCDs[obsID]['MOS1'])
            mos2CCD = str(self.sourceCCDs[obsID]['MOS2'])
            
            subprocess.run("cd "+workdir+";"+ \
                                ". $HEADAS/headas-init.sh;"+ \
                                ". $SAS_DIR/setsas.sh;"+ \
                                '''export SAS_CCF="`pwd`/ccf.cif";'''+ \
                                ";", shell=True)
        return True


##-- the main function --##
def main (args):

    #-- create an object of class spectra --#
    obj = epicPileup(ra=args.ra, dec=args.dec, workdir=args.workdir, \
                     sas_dir=args.sas_dir, headas=args.headas, sas_ccfpath=args.sas_ccfpath)
    
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

    #-- run the spectra functions --#
    obj.readPickleFile()
    obj.checkPileUp()

    print('\nHurray! The Spectra method ran successfully.')
    if len(obj.smallMode) != 0:
        print(f'\nThe following obsIDs have Small-mode MOS data.\n{obj.smallMode}')
    if len(obj.badObs) != 0:
        print('These obsIDs were excluded from analysis: ', obj.badObs)

    return True



if __name__=="__main__":
    
    description = 'Program to check and correct for EPIC Pile-up issue.'
    
    parser = argparse.ArgumentParser(description=description)   #--create a ArgumentParser object.

    #-- general arguments --#
    parser.add_argument('--ra', action='store', type=float, default=342.567, \
                        help='right ascension of the object. (default:%(default)s, AT-2018fyk)')
    parser.add_argument('--dec', action='store', type=float, default=-44.86, \
                        help='declination of the object. (default:%(default)s, AT-2018fyk)')   
    parser.add_argument('--workdir', action='store', type=str, default='/media/suyog/DATA/xmm_obs', \
                        help='directory where obsid folders will be stored. (default:%(default)s)')
    parser.add_argument('--obsIDs', nargs='+', action='store', default=None, #['0831790201'], \
                        help='''obsIDs for which some specific function has to executed. 
                                Valid only when --method argument is specified. (default:%(default)s)''')
    parser.add_argument('--objName', action='store', default=None, \
                        help='name of the obj used to locate the xmmObj pickle file.')
    
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

    #-- parse the arguments --#
    args = parser.parse_args()   #--parse all the arguments.
    #print(args)

    #-- call the main function --#
    main(args)



#################### End of Program #########################
#############################################################




