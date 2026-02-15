
####################################
###      EPIC Pile-up Issue      ###
####################################


import os
import subprocess
import logging

import glob
import argparse

from epicObj import epicObj

##-- the EPIC Pile-up class --##
class epicPileup (epicObj):

    ##-- check for Pile-up --##
    def checkPileUp (self):
        """
        Uses `epatplot` to check for the Pile-up issue.

        First a GTI filtered Source Event list is extracted out of `PN_CCD#.evts`
        and then `epatplot` is run to check the observed Spectra with the expected
        pattern distribution function.

        `epatplot` calculates two diagnostic numbers which may be used to assess 
        the presence of pile-up: In the absence of pile-up, the 0.5 - 2.0 keV (default) 
        observed-to-model singles and doubles pattern fractions ratios should both be 
        consistent with 1.0 within statistical errors (1 sigma errors are given). 
        If pile-up is present, the singles ratio will be smaller than 1.0 and the 
        doubles ratio will be larger than 1.0.
        See: `https://www.cosmos.esa.int/web/xmm-newton/sas-thread-epatplot`

        Note that `pnGTI.fits` GTI file is used instead of `combinedGTI_obsID.fits`
        because the latter results in very few useful data points.

        If the obsID is found to be Piled-up then the user should run the `correctPileUp`
        function to correct for it.

        :Input: Concatenated and Calibrated Event list, pnGTI file, 
               `sourceCCDs` and the Source location parameters.
        :Output: The Source Annulus radii together with Filtered Event List
                `source_PN_filtered_obsID.evts` and the corresponding image 
                `source_PN_filteredPattern.ps`.

        :Note: The term `filtered` means that the GTI have been applied.
        """
        logging.info('\nStarting Pile-up Checking.')

        #-- set the environment variables --#
        os.environ['SAS_DIR'] = self.sas_dir
        os.environ['HEADAS'] = self.headas
        os.environ['SAS_CCFPATH'] = self.sas_ccfpath
        
        #-- iterating for all the obsIDs --#      
        for obsID in self.obsIDs:
            logging.info(f'\nChecking Pile-up in {obsID}.')
            workdir = self.workdir+'/'+obsID+'/work'
            
            #-- get the location parameters --#
            locParams = self.sourceLoc[obsID]
            srcX, srcY, srcR = str(locParams['x']), str(locParams['y']), str(locParams['r'])
            #print(srcX, srcY, srcR)

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
                               " plotfile=source_PN_filteredPattern.ps;"
                           , shell=True)
            if self.showFig:
                subprocess.run("cd "+workdir+";"+ \
                            "evince source_PN_filteredPattern.ps;"
                            , shell=True)
            logging.info(f'\nChecked Pile-up for {obsID}.')

        return logging.info('\nPile-up checking finished.')


    ##-- correct the Pile-up --##
    def correctPileUp (self, srcRin=None):
        """
        Uses `epatplot` to correct and check for the Pile-up issue by removing
        some pixels from the center of the Source circle.

        The observation is corrected for Pile-up by removing some pixels from 
        the center of the Source. Now, this checking for Pile-up need be done
        only once, for PN images, since, it will be the same for MOS, probably.

        So, if Pile-up is found, one has to manually check for by iteratively 
        removing more and more center pixels and find when the Pile-up becomes 
        negligible. Then the Source circle would be an Annulus for this region.
        The same Annulus parameters can be repeated for MOS1&2.
        See: `https://www.cosmos.esa.int/web/xmm-newton/sas-thread-epatplot`

        Note that `pnGTI.fits` GTI file is used instead of `combinedGTI_obsID.fits`
        because the latter results in very few useful data points.

        The result of this procedure would be the radii values of the Source
        Annulus for which Pile-up is negligible. These new found `rIn` and `rOut`
        are then saved into to the `epicObj` and are written to `ccd_coords_info.pickle`
        file. Henceforth, for these Pile-up obsIDs the Source region defined
        from `rIn` to `rOut` will be used for obtaing the Spectra.

        :Input: Concatenated and Calibrated Event list, pnGTI file, 
               `sourceCCDs` and the Source location parameters.
        :Output: The Source Annulus radii together with Filtered Event List
                `source_PN_filtered_obsID.evts` and the corresponding image 
                `source_PN_filteredPattern.ps`, and the Annulus files 
                `source_PN_filteredAnnulus_obsID.evts` and `source_PN_filteredPattern_Annulus.ps` 
                showing resolution of the Pile-up issue.

        :Note: The term `filtered` means that the GTI have been applied.
        """
        logging.info('\nStarting Pile-up Correction.')
        if srcRin is None:
            return logging.error('Please provide a valid inner radius.')
        elif srcRin == 0:      
            for obsID in self.obsIDs:
                self.sourceLoc[obsID]['rIn'] = None
                logging.info('\nInner radius of the Source circle has been changed to None.')
            return logging.info(self.sourceLoc[obsID])
        else:
            srcRin = str(srcRin)

        #-- set the environment variables --#
        os.environ['SAS_DIR'] = self.sas_dir
        os.environ['HEADAS'] = self.headas
        os.environ['SAS_CCFPATH'] = self.sas_ccfpath
        
        #-- iterating for all the obsIDs --#      
        for obsID in self.obsIDs:
            logging.info('\nCorrecting Pile-up in {}.'.format(obsID))
            workdir = self.workdir+'/'+obsID+'/work'
            
            #-- get the location parameters --#
            locParams = self.sourceLoc[obsID]
            srcX, srcY, srcR = str(locParams['x']), str(locParams['y']), str(locParams['r'])
            logging.info(f'Source X, Y and R: {srcX} {srcY} {srcR}')

            #-- get the CCD numbers --#
            pnCCD = str(self.sourceCCDs[obsID]['PN'])
            mos1CCD = str(self.sourceCCDs[obsID]['MOS1'])
            mos2CCD = str(self.sourceCCDs[obsID]['MOS2'])
            
            #-- correct the Pile-up in PN --#
            subprocess.run("cd "+workdir+";"+ \
                           ". $HEADAS/headas-init.sh;"+ \
                           ". $SAS_DIR/setsas.sh;"+ \
                           '''export SAS_CCF="`pwd`/ccf.cif";'''+ \
                           "evselect table=PN_CCD"+pnCCD+".evts"+ \
                               " withfilteredset=yes filteredset=source_PN_filteredAnnulus_"+obsID+".evts"+ \
                               " keepfilteroutput=yes"+ \
                               " expression='((X,Y) in ANNULUS("+srcX+","+srcY+","+srcRin+","+srcR+"))"+ \
                               " && gti(pnGTI.fits,TIME)';"+ \
                           "epatplot set=source_PN_filteredAnnulus_"+obsID+".evts"+ \
                               " plotfile=source_PN_filteredAnnulus_Pattern.ps;"
                           , shell=True)
            if self.showFig:
                subprocess.run("cd "+workdir+";"+ \
                            "evince source_PN_filteredAnnulus_Pattern.ps;"
                            , shell=True)

            #-- save the new radius parameters --#
            self.sourceLoc[obsID]['rIn'] = srcRin
            self.sourceLoc[obsID]['rOut'] = srcR
            logging.info(f'\nCorrected Pile-up for {obsID}.')

        return logging.info('\nPile-up correction finished.')


    ##-- save the results --##
    def save_pileUpResults (self):
        """
        Saves the `*filteredPattern.ps` and `*filteredAnnulus*.ps`
        to the results directory.
        """
        logging.info('\nSaving the Pile-up results for each obsID to a common results folder.')

        maindir = self.workdir
        #objName+".dat;"
        #-- make the results directory --#
        resultdir = maindir+'/results/epic-Pileup'
        if not os.path.isdir(resultdir):
            if not os.path.isdir(maindir+'/results'):
                subprocess.run("cd "+maindir+";"+ \
                               "mkdir results/", shell=True)
            subprocess.run("cd "+maindir+"/results;"
                           "mkdir "+resultdir+";", shell=True)

        #-- iterating for all the obsIDs --#      
        for obsID in self.obsIDs:
            logging.info(f'\nSaving results for obsID {obsID}.')
            workdir = maindir+'/'+obsID+'/work'
            
            #-- copy results --#
            file1 = glob.glob(workdir+'/*filteredAnnulus*.ps')
            file2 = glob.glob(workdir+'/*filteredPattern.ps')
            files = [file[0] for file in [file1, file2] if len(file) != 0]

            for file in files:
                fname = os.path.basename(file)         #--get file name from the glob path.
                fname = fname.replace('_'+obsID, '')   #--remove obsID from file name, if it is already there.
                fname = obsID +'_'+ fname              #--add the obsID at the start of file name.
                
                subprocess.run("cp "+file+" "+resultdir+"/"+fname, shell=True)

            logging.info(f'Save for obsID {obsID} done.')

        return logging.info('\nSaved the result!')


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
            logging.error('\nNo obsIDs found for the given location. Confirm that you are connected to the Internet!')
    else:
        obj.objName = args.objName

    #-- check if obsID has been passed --#
    if args.obsIDs is not None:
        obj.obsIDs = args.obsIDs
        logging.info('Using the obsIDs passed.')

    if args.saveResults:
        obj.save_pileUpResults()
        return True

    obj.showFig = args.showFig

    #-- run the functions --#
    obj.readCCDcoordsPickle()
    if args.correctPileUp:
        obj.correctPileUp(args.srcRin)
        obj.writeCCDcoordsPickle()
    else:
        obj.checkPileUp()
    obj.save_pileUpResults()

    logging.info('\nHurray! Pile-up correction ran successfully.')
    if len(obj.smallMode) != 0 and args.obsIDs is None:
        logging.info(f'\nThe following obsIDs have Small-mode MOS data.\n{obj.smallMode}')
    if len(obj.badObs) != 0:
        logging.info(f'These obsIDs were excluded from analysis: {obj.badObs}')

    return True


##-------------------------------------------------------------------------------------------##

if __name__=="__main__":
    
    description = 'Program to check and correct for EPIC Pile-up issue.'
    
    parser = argparse.ArgumentParser(description=description)   #--create a ArgumentParser object.

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

    #-- flags for `correctPileUp` function --#
    parser.add_argument('--correctPileUp', action='store_true', default=False, \
                        help='run correctPileUp function. (default:%(default)s)')
    parser.add_argument('--srcRin', action='store', type=int, default=None, \
                        help='the inner radius of the Source circle, in Sky Coordinates. \
                             (default:%(default)s)')
    
    parser.add_argument('--showFig', action='store_true', default=False, \
                        help='show the epatplot output in evince. (default:%(default)s)')
    parser.add_argument('--saveResults', action='store_true', default=False, \
                        help='run only save_pileUpResults function. (default:%(default)s)')

    #-- parse the arguments --#
    args = parser.parse_args()   
    #print(args)

    #-- call the main function --#
    main(args)



#################### End of Program #########################
#############################################################




