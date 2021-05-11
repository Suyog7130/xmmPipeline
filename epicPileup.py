
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

This procedure will give the Source Annulus radii for which Pile-up is 
negligible for any particular obsIDs. Henceforth, it is this Source region
for which the Light Curve and the Spectra extraction should be done.

So, I will have to: 
    -> Modify `xmmPipeline.py` so that it knows which obsIDs have the 
       Pile-up issue and need the Source region to lie within an Annulus.
    -> Modify `spectra.py` so that the Source Spectra are extracted from 
       within this Annulus.

See the notes of Google Doc for open questions.

09 May 2021:
---
`epatplot` calculates two diagnostic numbers which may be used to assess 
the presence of pile-up: In the absence of pile-up, the 0.5 - 2.0 keV (default range) 
observed-to-model singles and doubles pattern fractions ratios should both be 
consistent with 1.0 within statistical errors (1 sigma errors are given). 
If pile-up is present, the singles ratio will be smaller than 1.0 and the 
doubles ratio will be larger than 1.0.
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


##-- the EPIC Pile-up class --##
class epicPileup (epicObj):

    ##-- check for Pile-up --##
    def checkPileUp (self):
        """
        Uses `epatplot` to check and correct for the Pile-up issue by removing
        some pixels from the center of the Source circle.

        The observation is corrected for Pile-up by removing some pixels from 
        the center of the Source. Now, this checking for Pile-up need be done
        only once, for PN images, since, it will be the same for MOS, probably.

        So, if Pile-up is found, one has to manually check for by iteratively 
        removing more and more center pixels and find when the Pile-up becomes 
        negligible. Then the Source circle would be an Annulus for this region.
        The same Annulus parameters can be repeated for MOS1&2.

        Note that `pnGTI.fits` GTI file is used instead of `combinedGTI_obsID.fits`
        because the latter results in very few useful data points.

        The result of this procedure would be the radii values of the Source
        Annulus for which Pile-up is negligible. Henceforth, this Source region
        will have to be used for obtaing the Spectra.

        Input: Concatenated and Calibrated Event list, pnGTI file, 
               `sourceCCDs` and the Source location parameters.
        Output: The Source Annulus radii together with Filtered Event List
                `source_PN_filtered_obsID.evts` and the corresponding image 
                `source_PN_filteredPattern.ps`, and the Annulus files 
                `source_PN_filteredAnnulus_obsID.evts` and `source_PN_filteredPattern_Annulus.ps` 
                showing resolution of the Pile-up issue.

        Note: The term `filtered` means that the GTI have been applied.
        """

        #-- set the environment variables --#
        os.environ['SAS_DIR'] = self.sas_dir
        os.environ['HEADAS'] = self.headas
        os.environ['SAS_CCFPATH'] = self.sas_ccfpath
        
        #-- iterating for all the obsIDs --#      
        for obsID in self.obsIDs:
            print('\nChecking Pile-up in {}.'.format(obsID))
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
            
            #-- check for Pile-up in PN --#
            subprocess.run("cd "+workdir+";"+ \
                           ". $HEADAS/headas-init.sh;"+ \
                           ". $SAS_DIR/setsas.sh;"+ \
                           '''export SAS_CCF="`pwd`/ccf.cif";'''+ \
                           "evselect table=PN_CCD"+pnCCD+".evts"+ \
                               " withfilteredset=yes filteredset=source_PN_filtered_"+obsID+".evts"+ \
                               " keepfilteroutput=yes"+ \
                               " expression='((X,Y) in CIRCLE("+srcX+","+srcY+","+srcR+"))"+ \
                               " && gti(pnGTI.fits,TIME)';"+ \
                           "epatplot set=source_PN_filtered_"+obsID+".evts"+ \
                               " plotfile=source_PN_filteredPattern.ps;"+ \
                           "evince source_PN_filteredPattern.ps;"
                           , shell=True)
            print(f'\nCorrected Pile-up for {obsID}.')

        return print('\nPile-up correction finished.')


##-------------------------------------------------------------------------------------------##

##-- the main function --##
def main (args):

    #-- create an object of epicPileup class --#
    obj = epicPileup(ra=args.ra, dec=args.dec, workdir=args.workdir, \
                     sas_dir=args.sas_dir, headas=args.headas, sas_ccfpath=args.sas_ccfpath)
    
    #-- get obsIDs and the objName --#
    if args.objName is None:
        try:
            obj.findObsIDs()
        except KeyError:
            print('\nNo obsIDs found for the given location. Confirm that you are connected to the Internet!')
    else:
        obj.objName = args.objName

    #-- check if obsID has been passed --#
    if args.obsIDs is not None:
        obj.obsIDs = args.obsIDs
        print('Using the obsIDs passed.')

    #-- run the functions --#
    obj.readCCDcoordsPickle()
    obj.checkPileUp()
    #obj.correctPileUp()

    print('\nHurray! Pile-up correction ran successfully.')
    if len(obj.smallMode) != 0 and args.obsIDs is None:
        print(f'\nThe following obsIDs have Small-mode MOS data.\n{obj.smallMode}')
    if len(obj.badObs) != 0:
        print('These obsIDs were excluded from analysis: ', obj.badObs)

    return True


##-------------------------------------------------------------------------------------------##

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





