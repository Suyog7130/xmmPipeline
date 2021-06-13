
#######################################
###      Spectra of XMM Object      ###
#######################################

"""
10th May 2021:
---
Am making several changes to the code procedure arrangements. 
See the GitHub repo, the Notes on Google Docs and the documentation for more information.
For previous docstring comments, see earlier code files, namely `spectra.py`

11th May 2021:
---
Since model fitting for individual obsIDs has to be done manually, I don't think
`xspec_fitSpectra` would be used now. `pyXspec.py` will have to be run manually
for each obsIDs with the different model parameters varrying.

20th May 2021:
---
Completing the pending work of fitting Models to the Spectra.
Two criterions to be checked for each obsIDs:
    - whether obsID is Piled-up? DONE!
    - whether obsID is in Small-mode? DONE!

21st May 2021:
---
Changed to using Flare Background filtered Event Lists for Spectra extraction.
Found a way to save the output parameter values from `pyXspec.py`.
Completing the SpecModel fitting now.

13th June 2021:
---
Completing the work for obtaining individual instrument bkgCircs.
    * Adding `runIndiBkgCircFuncs` to the main function.
"""

import os
import subprocess
import requests
import wget

import json
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



##-- the epicSpectra object class --##
class epicSpectra (epicObj):

    def __init__ (self, doNotOverwriteModel=False, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.doNotOverwriteModel = doNotOverwriteModel


    ##-- update location parameters in the pickle file --##
    def updateCCDcoordsPickle (self):
        """
        Writes the background location parameter dictionary for individual instruments
        ``backgroundLoc_indi``  to `ccd_coords_info.pickle` file.
        If the file is already present in the directory, then the values for the 
        corresponding keys are updated.

        :Input: ``backgroundLoc_indi`` dictionary containing the location parameters.
        :Output: Updated ``ccd_coords_info.pickle`` file.

        NOTES
        -----
        If no CCDcoordsPickle file is available in the ``workdir``, then the whole
        file is created with other values been saved as well.
        """
        maindir = self.workdir
        print('\nWriting indi inst bkgCircs data to ccd_coords_info.pickle file.\n')

        #-- iterating for all the obsIDs --#      
        for obsID in self.obsIDs:
            print('Saving pickle for obsID {}.'.format(obsID))
            workdir = maindir+'/'+obsID+'/work'
            fname = workdir+'/'+'ccd_coords_info.pickle'

            #-- load the `ccd_coords_info.pickle` if already present --#
            if os.path.isfile(fname):
                result = pickle.load(open(fname, 'rb'))
                    #-- create a backup file --#
                subprocess.run("cd "+workdir+";"+ \
                               "cp ccd_coords_info.pickle ccd_coords_info.bak;"
                               , shell=True)
            else:
                result = {}

            #-- save backgroundLoc_indi dictionary --#
            dic = result.get('backgroundLoc_indi', None)
            if dic is None:
                result['backgroundLoc_indi'] = {}
                for inst in 

            
            #-- save the CCD and coords info in a pickle file --#
            result['sourceCCDs'] = self.sourceCCDs[obsID]
            result['sourceLoc'] = self.sourceLoc[obsID]
            result['backgroundLoc'] = self.backgroundLoc[obsID]
            result['otherSources'] = self.otherSources[obsID]
            result['smallMode'] = self.smallMode[obsID]
            if obsID in self.badObs:
                result['badObs'] = True
            else:
                result['badObs'] = False

            outfile = open(fname, 'wb')
            pickle.dump(result, outfile)
            outfile.close()

        return print('\nWrote location parameters to the pickle file.')


    ##-- extract the Spectra in Image Mode --##
    def extractSpectra_imageMode (self):
        """
        Extracts the MOS12 and PN spectra in Image Mode.
        See: https://www.cosmos.esa.int/web/xmm-newton/sas-thread-mos-spectrum,
             https://www.cosmos.esa.int/web/xmm-newton/sas-thread-pn-spectrum

        :Input: Flared Background filtered PN, MOS1&2 Event Lists, `PNclean.ds`, `MOS1clean.ds` and
               `MOS2clean.ds`. Alongwith SrcBkg coordinates and SourceCCDs resulting from `epicPipeline.py`
        :Output: MOS12 and PN SrcBkg Spectra in Image Mode, a Redistribution Matrix (`rmf`) file and
                an Effective Area Vector (`arf`) file.

        Note that for Small Mode obsIDs, since the background is taken from PN image data,
        instead of the overlap region data b'cuz MOS12 region was of small size, the background
        region will lie outside of the image. Thus, only the Source Spectra can possibly be 
        extracted for Small Mode obsIDs. However, the 'specgroup' command requires both the 
        Source and the Background Spectra in order to create the grouped Spectrum. Therefore,
        as of now, I am completely skipping getting the MOS12 Spectra for Small mode obsIDs.

        20th May 2021: Now works with Pile-up corrected obsIDs.
        Separate Source Spectrum file is generated for Piled-up cases, `inst_spectrum_source_annulus_obsID.fits`
        All the other file names, including the Grouped Spectra is not altered.
        """
        print('\nStarting Spectra extraction for MOS12 and PN in Image Mode.')
        
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
            
            #-- Piled-up cases --#
            srcRin = locParams.get('rIn', None)
            if not self.ignorePileup and srcRin is not None:
                srcRin = str(srcRin)
                srcRout = str(locParams['rOut'])
                print('\nThis observation was piled-up. Using ANNULUS for Source.')
                srcSpectrumSet = "spectrum_source_annulus_"+obsID+".fits"
                srcFilterExp = "((X,Y) in ANNULUS("+srcX+","+srcY+","+srcRin+","+srcRout+"))';"
            else:
                srcSpectrumSet = "spectrum_source_"+obsID+".fits"
                srcFilterExp = "((X, Y) IN circle("+srcX+","+srcY+","+srcR+"))';"

            #-- extract MOS1 Spectra --#
            if self.smallMode[obsID]:
                print(f'\n{obsID} is in Small Mode. Skipping MOS1 Spectra extraction for it.')
            else:
                subprocess.run("cd "+workdir+";"+ \
                               ". $HEADAS/headas-init.sh;"+ \
                               ". $SAS_DIR/setsas.sh;"+ \
                               '''export SAS_CCF="`pwd`/ccf.cif";'''+ \
                               "evselect table=MOS1clean.ds"+ \
                                   " withspectrumset=yes spectrumset=MOS1_"+srcSpectrumSet+ \
                                   " energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999"+ \
                                   " expression='#XMMEA_EM && (PATTERN<=12) && "+srcFilterExp+ \
                               "evselect table=MOS1clean.ds"+ \
                                   " withspectrumset=yes spectrumset=MOS1_spectrum_background_"+obsID+".fits"+ \
                                   " energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999"+ \
                                   " expression='#XMMEA_EM && (PATTERN<=12) && ((X,Y) in CIRCLE("+ \
                                   Bx1+","+By1+","+Br1+"))||((X,Y) in CIRCLE("+Bx2+","+By2+","+Br2+"))';"+ \
                               "backscale spectrumset=MOS1_"+srcSpectrumSet+" badpixlocation=MOS1clean.ds;"+ \
                               "backscale spectrumset=MOS1_spectrum_background_"+obsID+".fits badpixlocation=MOS1clean.ds;"+ \
                               "rmfgen spectrumset=MOS1_"+srcSpectrumSet+" rmfset=MOS1_"+obsID+".rmf;"+ \
                               "arfgen spectrumset=MOS1_"+srcSpectrumSet+ \
                                   " arfset=MOS1_"+obsID+".arf withrmfset=yes rmfset=MOS1_"+obsID+".rmf"+ \
                                   " badpixlocation=MOS1clean.ds detmaptype=psf;"+ \
                               "specgroup spectrumset=MOS1_"+srcSpectrumSet+ \
                                   " mincounts=20 oversample=3 rmfset=MOS1_"+obsID+".rmf"+ \
                                   " arfset=MOS1_"+obsID+".arf backgndset=MOS1_spectrum_background_"+obsID+".fits"+ \
                                   " groupedset=MOS1_spectrum_grouped_"+obsID+".fits;" 
                               #"fv MOS1_spectrum_grouped_"+obsID+".fits;"
                               , shell=True)
                print('\nMOS1 Spectra extracted.')

            #-- extract MOS2 Spectra --#
            if self.smallMode[obsID]:
                print(f'\n{obsID} is in Small Mode. Skipping MOS2 Spectra extraction for it.')
            else:
                subprocess.run("cd "+workdir+";"+ \
                               ". $HEADAS/headas-init.sh;"+ \
                               ". $SAS_DIR/setsas.sh;"+ \
                               '''export SAS_CCF="`pwd`/ccf.cif";'''+ \
                               "evselect table=MOS2clean.ds"+ \
                                   " withspectrumset=yes spectrumset=MOS2_"+srcSpectrumSet+ \
                                   " energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999"+ \
                                   " expression='#XMMEA_EM && (PATTERN<=12) && "+srcFilterExp+ \
                               "evselect table=MOS2clean.ds"+ \
                                   " withspectrumset=yes spectrumset=MOS2_spectrum_background_"+obsID+".fits"+ \
                                   " energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999"+ \
                                   " expression='#XMMEA_EM && (PATTERN<=12) && ((X,Y) in CIRCLE("+ \
                                   Bx1+","+By1+","+Br1+"))||((X,Y) in CIRCLE("+Bx2+","+By2+","+Br2+"))';"+ \
                               "backscale spectrumset=MOS2_"+srcSpectrumSet+" badpixlocation=MOS2clean.ds;"+ \
                               "backscale spectrumset=MOS2_spectrum_background_"+obsID+".fits badpixlocation=MOS2clean.ds;"+ \
                               "rmfgen spectrumset=MOS2_"+srcSpectrumSet+" rmfset=MOS2_"+obsID+".rmf;"+ \
                               "arfgen spectrumset=MOS2_"+srcSpectrumSet+ \
                                   " arfset=MOS2_"+obsID+".arf withrmfset=yes rmfset=MOS2_"+obsID+".rmf"+ \
                                   " badpixlocation=MOS2clean.ds detmaptype=psf;"+ \
                               "specgroup spectrumset=MOS2_"+srcSpectrumSet+ \
                                   " mincounts=20 oversample=3 rmfset=MOS2_"+obsID+".rmf"+ \
                                   " arfset=MOS2_"+obsID+".arf backgndset=MOS2_spectrum_background_"+obsID+".fits"+ \
                                   " groupedset=MOS2_spectrum_grouped_"+obsID+".fits;" 
                               #"fv MOS2_spectrum_grouped_"+obsID+".fits;"
                               , shell=True)
                print('\nMOS2 Spectra extracted.')

            #-- extract PN Spectra --#
            subprocess.run("cd "+workdir+";"+ \
                           ". $HEADAS/headas-init.sh;"+ \
                           ". $SAS_DIR/setsas.sh;"+ \
                           '''export SAS_CCF="`pwd`/ccf.cif";'''+ \
                           "evselect table=PNclean.ds"+ \
                               " withspectrumset=yes spectrumset=PN_"+srcSpectrumSet+ \
                               " energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479"+ \
                               " expression='#XMMEA_EP && (FLAG==0) && (PATTERN<=4) && "+srcFilterExp+ \
                           "evselect table=PNclean.ds"+ \
                               " withspectrumset=yes spectrumset=PN_spectrum_background_"+obsID+".fits"+ \
                               " energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479"+ \
                               " expression='#XMMEA_EP && (FLAG==0) && (PATTERN<=4) && ((X,Y) in CIRCLE("+ \
                               Bx1+","+By1+","+Br1+"))||((X,Y) in CIRCLE("+Bx2+","+By2+","+Br2+"))';"+ \
                           "backscale spectrumset=PN_"+srcSpectrumSet+" badpixlocation=PNclean.ds;"+ \
                           "backscale spectrumset=PN_spectrum_background_"+obsID+".fits badpixlocation=PNclean.ds;"+ \
                           "rmfgen spectrumset=PN_"+srcSpectrumSet+" rmfset=PN_"+obsID+".rmf;"+ \
                           "arfgen spectrumset=PN_"+srcSpectrumSet+ \
                               " arfset=PN_"+obsID+".arf withrmfset=yes rmfset=PN_"+obsID+".rmf"+ \
                               " badpixlocation=PNclean.ds detmaptype=psf;"+ \
                           "specgroup spectrumset=PN_"+srcSpectrumSet+ \
                               " mincounts=20 oversample=3 rmfset=PN_"+obsID+".rmf"+ \
                               " arfset=PN_"+obsID+".arf backgndset=PN_spectrum_background_"+obsID+".fits"+ \
                               " groupedset=PN_spectrum_grouped_"+obsID+".fits;"
                           #"fv PN_spectrum_grouped_"+obsID+".fits;"
                           , shell=True)
            print('\nPN Spectra extracted.')
            print('\nExtracting Image Mode Spectra for obsID {} finished.'.format(obsID)) 

        return print('\nCompleted extracting Spectra in Image Mode.')


    ##-- fit Spectra to extracted group Spectra data --##
    def xspec_fitSpectra (self, instName='PN'):
        """
        Fits a model Spectra to the extracted Spectra data using `pyXspec` package.
        Since, `pyXspec` requires HEA initialisation, a separate routine `pyXspec.py`
        is required to run it.
        See: https://www.cosmos.esa.int/web/xmm-newton/sas-thread-xspec

        :Input: Extracted Spectra FITS file and input model parameters.
        :Output: Spectra plots and `specModelParams.json` file.

        20th May 2021: Adding feature to allow model parameter value inputs.
        """
        print('\nStarting to fit a model Spectra to the extracted Spectra using xspec.')
        
        #-- set the environment variables --#
        os.environ['SAS_DIR'] = self.sas_dir
        os.environ['HEADAS'] = self.headas
        os.environ['SAS_CCFPATH'] = self.sas_ccfpath
        
        #-- iterating for all the obsIDs --#      
        for obsID in self.obsIDs:
            print('\nFitting Spectra for obsID {}.'.format(obsID))
            workdir = self.workdir+'/'+obsID+'/work'

            smallMode = self.smallMode[obsID]
            smallMode = strToBool(smallMode, inverse=True)
            showFig = strToBool(self.showFig, inverse=True)
            doNotOverwriteModel = strToBool(self.doNotOverwriteModel, inverse=True)

            #-- read a modelParams dict --#
            model = self.model.replace(')', '\)').replace('(', '\(')
            modelParams = "'" + str(self.modelParams).replace("'", '"') + "'"
            #print(model, modelParams)

            #-- run `pyXspec.py` routine --#
            subprocess.run("cd "+workdir+";"+ \
                            ". $HEADAS/headas-init.sh;"+ \
                            ". $SAS_DIR/setsas.sh;"+ \
                            #'''export SAS_CCF="`pwd`/ccf.cif";'''+ \
                            "python3 ~/Dropbox/Dheeraj@MIT_2020-21/pyXspec.py --obsID "+obsID+ \
                                " --instName all --showFig "+showFig+" --smallMode "+smallMode+ \
                                " --model "+model+" --modelParams "+modelParams+ \
                                " --doNotOverwriteModel "+doNotOverwriteModel+";"
                            , shell=True)
            print('\nSpectra fitting for obsID {} finished.'.format(obsID))
            
        return print('Fitted model Spectra to the extracted Spectra.')


    ##-- function to save the final results --##
    def save_spectraResults (self):
        """
        Does one tasks
                Copies the Spectra images from each obsID directory to a results 
                folder in the main directory.
        """
        print('\nLastly saving results for each obsID to a common results folder.')

        maindir = self.workdir
        #objName+".dat;"
        #-- make the results directory --#
        resultdir = maindir+'/results/epic-Spectra'
        if not os.path.isdir(resultdir):
            if not os.path.isdir(maindir+'/results'):
                subprocess.run("cd "+maindir+";"+ \
                               "mkdir results/", shell=True)
            subprocess.run("cd "+maindir+"/results;"
                           "mkdir "+resultdir+";", shell=True)

        #-- iterating for all the obsIDs --#      
        for obsID in self.obsIDs:
            print('\nSaving results for obsID {}.'.format(obsID))
            workdir = maindir+'/'+obsID+'/work'
            
            #-- copy Event lists and other results --#
            spectraImgFiles = glob.glob(workdir+'/*spectra*.png')
            srcSpectra = glob.glob(workdir+'/*spectrum_source*.fits')
            bkgSpectra = glob.glob(workdir+'/*spectrum_background*.fits')
            rmfFiles = glob.glob(workdir+'/*.rmf')
            arfFiles = glob.glob(workdir+'/*.arf')
            groupedSpectra = glob.glob(workdir+'/*spectrum_grouped*.fits')
            paramsFile = glob.glob(workdir+'/*specModelParams.json')

            spectraFiles = spectraImgFiles + srcSpectra + bkgSpectra + rmfFiles + \
                           arfFiles + groupedSpectra + paramsFile

            for file in spectraFiles:
                fname = os.path.basename(file)         #--get file name from the glob path.
                fname = fname.replace('_'+obsID, '')   #--remove obsID from file name, if it is already there.
                fname = obsID +'_'+ fname              #--add the obsID at the start of file name.
                
                file = file.replace('(','\(').replace(')','\)')
                fname = fname.replace('(','\(').replace(')','\)')

                subprocess.run("cp "+file+" "+resultdir+"/"+fname, shell=True)

            print(f'Save for obsID {obsID} done.')

        return print('\nSaved the result!')


##-------------------------------------------------------------------------------------------##


##-- the main function --##
def main (args):
    
    #-- create an object of class spectra --#
    obj = epicSpectra(ra=args.ra, dec=args.dec, workdir=args.workdir, \
                      sas_dir=args.sas_dir, headas=args.headas, sas_ccfpath=args.sas_ccfpath, \
                      saveFig=args.noSaveFig, showFig=args.showFig, \
                      ignorePileup=args.ignorePileup, doNotOverwriteModel=args.doNotOverwriteModel)
    
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
    
    if args.saveResults:
        obj.save_spectraResults()
        return True

    #-- if ``runIndiBkgCircFuncs``, then run them and exit --#
    if args.inst is not None:
        obj.readCCDcoordsPickle()
        obj.otherSources_indi(inst=args.inst)
        obj.getBackgroundCircles_indi(inst=args.inst)
        obj.updateCCDcoordsPickle(inst=args.inst)
        return print('\nRan the IndiBkgCirc functions and obtained indi inst bkg circles. \
            The CCDcoordsPickle file has also been updated with backgroundLocIndi key.')

    #-- run the spectra functions --#
    obj.readCCDcoordsPickle()
    if args.method == 'extractSpectra':
        obj.extractSpectra_imageMode()
    if args.method == 'fitSpectra':
        obj.model, obj.modelParams = args.model, args.modelParams
        #print(obj.__dict__)
        obj.xspec_fitSpectra()
    obj.save_spectraResults()

    print('\nHurray! The Spectra method ran successfully.')
    #if len(obj.smallMode) != 0:
    #    print(f'\nThe following obsIDs have Small-mode MOS data.\n{obj.smallMode}')
    if len(obj.badObs)!=0:
        print('These obsIDs were excluded from analysis: ', obj.badObs)

    return True


if __name__=="__main__":
    
    description = 'Program to reduce RGS and OM data and obtain the Spectra of XMM objects.'
    
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
                        help='name of the obj used to locate the xmmObj pickle file. \
                              (default:%(default)s)')

    parser.add_argument('--inst', action='store', default=None, \
                        help='Spectral Analysis requires the three EPIC camera data \
                              to be separately processed. For this individual instrument \
                              background circles are required. Use this to pass the \
                              individual instrument for which bkgCircs have to be found. \
                              Also updates CCDcoordsPickle. (default:%(default)s)')
    
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

    #-- flags --#
    parser.add_argument('--noSaveFig', action='store_false', default=True, \
                        help='do not save the matplotlib plots? (default:%(default)s)')
    parser.add_argument('--showFig', action='store_true', default=False, \
                        help='show the matplotlib plots or not? (default:%(default)s)')
    parser.add_argument('--saveResults', action='store_true', default=False, \
                        help='run only save_spectraResults function. (default:%(default)s)')   
    parser.add_argument('--ignorePileup', action='store_true', default=False, \
                        help='ignore Pile-up in Piled-up obsIDs. (default:%(default)s)')
    parser.add_argument('--doNotOverwriteModel', action='store_true', default=False, \
                        help='save the specModelParams for the current model as a new \
                              dictionary in specModelParams.json. \
                              If the model is already present in the file, it is not \
                              overwritten by appending 1 to the new model name. \
                              (default:%(default)s)')

    #-- methods --#
    parser.add_argument('--method', action='store', type=str, default='extractSpectra', \
                        choices=['extractSpectra', 'fitSpectra'], \
                        help='what to do? `fitSpectra` calls pyXspec to fit models to \
                              the Spectra. (default:%(default)s)')

    #-- fitSpectra arguments --#
    parser.add_argument('--instName', action='store', default='PN', \
                        help='the instrument to use. (default:%(default)s)')
    parser.add_argument('--model', action='store', default='tbabs*zashift*(powerlaw)', \
                        help='name of the model to be used. (default:%(default)s)')
    parser.add_argument('--modelParams', action='store', type=eval, default={}, \
                        help='give a dictionary of model parameter values to use. \
                              Put double quotes for str values and enclose the dict \
                              within single quotes at the end.')

    #-- parse the arguments --#
    args = parser.parse_args()   #--parse all the arguments.
    #print('{}\n{}'.format(args, args.modelParams))
    args.runIndiBkgCircFuncs = strToBool(args.runIndiBkgCircFuncs)

    #-- call the main function --#
    main(args)
    





#################### End of Program #########################
#############################################################
