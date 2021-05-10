
#######################################
###      Spectra of XMM Object      ###
#######################################

"""
10th May 2021:
---
Am making several changes to the code procedure arrangements. 
See the GitHub repo, the Notes on Google Docs and the documentation for more information.
For previous docstring comments, see earlier code files, namely `spectra.py`
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


##-- the spectra object class --##
class spectra:

    ##-- initialize some common parameters --##
    def __init__(self, ra, dec, workdir, sas_dir, headas, sas_ccfpath):

        self.ra = ra
        self.dec = dec
        self.workdir = workdir
        
        self.sas_dir = sas_dir
        self.headas = headas
        self.sas_ccfpath = sas_ccfpath
        
        self.objName = None       #--name of the object found at the given location.
        self.obsIDs = list()
        self.badObs = list()      #--list of obsIDs, excluded from analysis, having total GTI below gti_combThreshold.
        self.sourceCCDs = {}      #--CCD numbers for each obsID.
        self.sourceLoc = {}       #--source location parameters for each obsID.
        self.backgroundLoc = {}   #--background circle coordinates and radius for each obsID.

        self.smallMode = {}       #--is obsID in smallMode, dict for all obsIDs.


    ##-- function to find the obsIDs --##
    def findObsIDs (self):
        """
        Input: Coordinates of the object and path to work directory.
        Output: List of obsIDs.
        """
        ra, dec, workdir = str(self.ra), str(self.dec), self.workdir
        print('\nLooking for obsIDs at RA={} and DEC={}\nWORKDIR is set at {}'.format(ra,dec,workdir))
        
        #-- check if workdir exists --#
        if not os.path.isdir(workdir):
            subprocess.run("sudo mkdir "+workdir, shell=True)
            
        #-- check if browse_extract_wget.pl file exists --#
        if not os.path.isfile(workdir+"/"+"browse_extract_wget.pl"):
            print('\nbrowse_extract_wget.pl not found.')
            subprocess.run("cd "+workdir+";"+
                           "sudo wget -q https://heasarc.gsfc.nasa.gov/FTP/heasarc/software/web_batch/browse_extract_wget.pl", shell=True)
            
            print('browse_extract_wget.pl downloaded.')
            print('\nPlease check the PERL path in the file. If required, correct the path given in first line and save the file.')
            subprocess.run("cd "+workdir+";"+
                           "sudo gedit browse_extract_wget.pl", shell=True)
            
        #-- download and save parts of xmmmaster table --##
        subprocess.run("cd "+workdir+";"+ \
                       "sudo chmod +x browse_extract_wget.pl;"+ \
                       "sudo ./browse_extract_wget.pl table=xmmmaster position='"+ ra+","+dec+ \
                       "' coordinates=equatorial outfile=obsIDs_list.dat", shell=True)
        
        #-- extract obsIDs from the file --#
        df = pd.read_csv(workdir+'/obsIDs_list.dat', sep='|', delim_whitespace=False, header=0)[:-1]  #--remove last line.
        cols = [s.strip() for s in df.columns.to_list()]  #--remove whitespace from column names.
        df.columns = cols

        obsIDs = []
        for i, s in enumerate(df['_Search_Offset'].fillna(0)):
            if s!=0 and float(s.split()[0].strip())<1.0:           #--remove those which are too far off.
                obsIDs.append( '0' + str(int(df['obsid'][i])) )    #--the database has a 0 at the start.
       
        self.obsIDs = obsIDs  

        #-- get the objName from the file --#
        objName = df.name.drop_duplicates().dropna().tolist()[0]
        self.objName = objName
        print(f'The object at (RA,DEC) = ({ra},{dec}) is {objName}')

        #-- copy and append objName to the file --#
        subprocess.run("cd "+workdir+";"+ \
                       "cp obsIDs_list.dat obsIDs_list_"+objName+".dat;", shell=True) 
        
        return print('\nFound '+str(len(obsIDs))+' obsIDs for the object at given position.')


    ##-- read the pickle file for location parameters --##
    def readPickleFile (self):
        """
        Reads the already obtained pickle file containing location parameters.
        """
        print('\nLoading the Pickle file obtained from xmmPipeline.')

        workdir, objName = self.workdir, self.objName
        fname = workdir+'/'+'output_xmmObj_'+objName+'.pickle'

        #-- check for the pickle file --#
        if not os.path.isfile(fname):
            return print(f'\nFile {fname} not found!')

        #-- load the epicObj --#
        epicObj = pickle.load(open(fname, 'rb'))

        #-- get the location and other parameters --#
        self.badObs = epicObj.badObs
        self.sourceCCDs = epicObj.sourceCCDs
        self.sourceLoc, self.backgroundLoc = epicObj.sourceLoc, epicObj.backgroundLoc
        self.smallMode, self.otherSources = epicObj.smallMode, epicObj.otherSources
        
        #-- remove badObs from obsIDs list to use --#
        obsIDs, badObs = set(self.obsIDs), set(self.badObs)
        self.obsIDs = list(obsIDs-badObs)

        return print('Read location parameters from the xmmPipeline Pickle file.')


    ##-- extract the Spectra in Image Mode --##
    def extractSpectra_imageMode (self):
        """
        Extracts the MOS12 and PN spectra in Image Mode.
        See: https://www.cosmos.esa.int/web/xmm-newton/sas-thread-mos-spectrum,
             https://www.cosmos.esa.int/web/xmm-newton/sas-thread-pn-spectrum

        Input: MOS12.evts, PN_CCD##.evts, SrcBkg coordinates and SourceCCDs obtained from xmmPipeline.py
        Output: MOS12 and PN SrcBkg Spectra in Image Mode, a Redistribution Matrix (rmf) file and
                an Effective Area Vector (arf) file.

        Note that for Small Mode obsIDs, since the background is taken from PN image data,
        instead of the overlap region data b'cuz MOS12 region was of small size, the background
        region will lie outside of the image. Thus, only the Source Spectra can possibly be 
        extracted for Small Mode obsIDs. However, the 'specgroup' command requires both the 
        Source and the Background Spectra in order to create the grouped Spectrum. Therefore,
        as of now, I am completely skipping getting the MOS12 Spectra for Small mode obsIDs.
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

            #-- get the CCD numbers --#
            pnCCD = str(self.sourceCCDs[obsID]['PN'])
            mos1CCD = str(self.sourceCCDs[obsID]['MOS1'])
            mos2CCD = str(self.sourceCCDs[obsID]['MOS2'])
            
            """
            #-- extract MOS12 Spectra --#
            if self.smallMode[obsID]:
                print(f'\n{obsID} is in Small Mode. Skipping MOS12 Spectra extraction for it.')
            else:
                subprocess.run("cd "+workdir+";"+ \
                               ". $HEADAS/headas-init.sh;"+ \
                               ". $SAS_DIR/setsas.sh;"+ \
                               '''export SAS_CCF="`pwd`/ccf.cif";'''+ \
                               #"ds9 MOS12_image.fits;"
                               "evselect table=MOS12.evts:EVENTS withspectrumset=yes spectrumset=MOS12_spectrum_source.fits"+ \
                                   " energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999"+ \
                                   " expression='#XMMEA_EM && (PATTERN<=12) && ((X, Y) IN circle("+ \
                                   srcX+","+srcY+","+srcR+"))';"+ \
                               "evselect table=MOS12.evts:EVENTS withspectrumset=yes spectrumset=MOS12_spectrum_background.fits"+ \
                                   " energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999"+ \
                                   " expression='#XMMEA_EM && (PATTERN<=12) && ((X,Y) in CIRCLE("+ \
                                   Bx1+","+By1+","+Br1+"))||((X,Y) in CIRCLE("+Bx2+","+By2+","+Br2+"))';"+ \
                               "backscale spectrumset=MOS12_spectrum_source.fits badpixlocation=MOS12.evts;"+ \
                               "backscale spectrumset=MOS12_spectrum_background.fits badpixlocation=MOS12.evts;"+ \
                               "rmfgen spectrumset=MOS12_spectrum_source.fits rmfset=MOS12.rmf;"+ \
                               "arfgen spectrumset=MOS12_spectrum_source.fits arfset=MOS12.arf withrmfset=yes rmfset=MOS12.rmf"+ \
                                   " badpixlocation=MOS12.evts detmaptype=psf;"+ \
                               "specgroup spectrumset=MOS12_spectrum_source.fits mincounts=25 oversample=3 rmfset=MOS12.rmf"+ \
                                   " arfset=MOS12.arf backgndset=MOS12_spectrum_background.fits"+ \
                                   " groupedset=MOS12_spectrum_grouped.fits;"
                               , shell=True)
                #self.xspec_fitSpectra(instName='MOS12')
                print('\nMOS12 Spectra extracted.')
            """

            #-- extract MOS1 Spectra --#
            if self.smallMode[obsID]:
                print(f'\n{obsID} is in Small Mode. Skipping MOS1 Spectra extraction for it.')
            else:
                subprocess.run("cd "+workdir+";"+ \
                               ". $HEADAS/headas-init.sh;"+ \
                               ". $SAS_DIR/setsas.sh;"+ \
                               '''export SAS_CCF="`pwd`/ccf.cif";'''+ \
                               "evselect table=MOS1_CCD"+mos1CCD+".evts:EVENTS"+ \
                                   " withspectrumset=yes spectrumset=MOS1_spectrum_source_"+obsID+".fits"+ \
                                   " energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999"+ \
                                   " expression='#XMMEA_EM && (PATTERN<=12) && ((X, Y) IN circle("+ \
                                   srcX+","+srcY+","+srcR+"))';"+ \
                               "evselect table=MOS1_CCD"+mos1CCD+".evts:EVENTS"+ \
                                   " withspectrumset=yes spectrumset=MOS1_spectrum_background_"+obsID+".fits"+ \
                                   " energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999"+ \
                                   " expression='#XMMEA_EM && (PATTERN<=12) && ((X,Y) in CIRCLE("+ \
                                   Bx1+","+By1+","+Br1+"))||((X,Y) in CIRCLE("+Bx2+","+By2+","+Br2+"))';"+ \
                               "backscale spectrumset=MOS1_spectrum_source_"+obsID+".fits badpixlocation=MOS1_CCD"+mos1CCD+".evts;"+ \
                               "backscale spectrumset=MOS1_spectrum_background_"+obsID+".fits badpixlocation=MOS1_CCD"+mos1CCD+".evts;"+ \
                               "rmfgen spectrumset=MOS1_spectrum_source_"+obsID+".fits rmfset=MOS1_"+obsID+".rmf;"+ \
                               "arfgen spectrumset=MOS1_spectrum_source_"+obsID+".fits"+ \
                                   " arfset=MOS1_"+obsID+".arf withrmfset=yes rmfset=MOS1_"+obsID+".rmf"+ \
                                   " badpixlocation=MOS1_CCD"+mos1CCD+".evts detmaptype=psf;"+ \
                               "specgroup spectrumset=MOS1_spectrum_source_"+obsID+".fits"+ \
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
                               "evselect table=MOS2_CCD"+mos2CCD+".evts:EVENTS"+ \
                                   " withspectrumset=yes spectrumset=MOS2_spectrum_source_"+obsID+".fits"+ \
                                   " energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999"+ \
                                   " expression='#XMMEA_EM && (PATTERN<=12) && ((X, Y) IN circle("+ \
                                   srcX+","+srcY+","+srcR+"))';"+ \
                               "evselect table=MOS2_CCD"+mos2CCD+".evts:EVENTS"+ \
                                   " withspectrumset=yes spectrumset=MOS2_spectrum_background_"+obsID+".fits"+ \
                                   " energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999"+ \
                                   " expression='#XMMEA_EM && (PATTERN<=12) && ((X,Y) in CIRCLE("+ \
                                   Bx1+","+By1+","+Br1+"))||((X,Y) in CIRCLE("+Bx2+","+By2+","+Br2+"))';"+ \
                               "backscale spectrumset=MOS2_spectrum_source_"+obsID+".fits badpixlocation=MOS2_CCD"+mos2CCD+".evts;"+ \
                               "backscale spectrumset=MOS2_spectrum_background_"+obsID+".fits badpixlocation=MOS2_CCD"+mos2CCD+".evts;"+ \
                               "rmfgen spectrumset=MOS2_spectrum_source_"+obsID+".fits rmfset=MOS2_"+obsID+".rmf;"+ \
                               "arfgen spectrumset=MOS2_spectrum_source_"+obsID+".fits"+ \
                                   " arfset=MOS2_"+obsID+".arf withrmfset=yes rmfset=MOS2_"+obsID+".rmf"+ \
                                   " badpixlocation=MOS2_CCD"+mos2CCD+".evts detmaptype=psf;"+ \
                               "specgroup spectrumset=MOS2_spectrum_source_"+obsID+".fits"+ \
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
                           #"ds9 PN_CCD"+pnCCD+"_image.fits -scale log;"
                           "evselect table=PN_CCD"+pnCCD+".evts:EVENTS"+ \
                               " withspectrumset=yes spectrumset=PN_spectrum_source_"+obsID+".fits"+ \
                               " energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479"+ \
                               " expression='#XMMEA_EP && (FLAG==0) && (PATTERN<=4) && ((X, Y) IN circle("+ \
                               srcX+","+srcY+","+srcR+"))';"+ \
                           "evselect table=PN_CCD"+pnCCD+".evts:EVENTS"+ \
                               " withspectrumset=yes spectrumset=PN_spectrum_background_"+obsID+".fits"+ \
                               " energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479"+ \
                               " expression='#XMMEA_EP && (FLAG==0) && (PATTERN<=4) && ((X,Y) in CIRCLE("+ \
                               Bx1+","+By1+","+Br1+"))||((X,Y) in CIRCLE("+Bx2+","+By2+","+Br2+"))';"+ \
                           "backscale spectrumset=PN_spectrum_source_"+obsID+".fits badpixlocation=PN_CCD"+pnCCD+".evts;"+ \
                           "backscale spectrumset=PN_spectrum_background_"+obsID+".fits badpixlocation=PN_CCD"+pnCCD+".evts;"+ \
                           "rmfgen spectrumset=PN_spectrum_source_"+obsID+".fits rmfset=PN_"+obsID+".rmf;"+ \
                           "arfgen spectrumset=PN_spectrum_source_"+obsID+".fits"+ \
                               " arfset=PN_"+obsID+".arf withrmfset=yes rmfset=PN_"+obsID+".rmf"+ \
                               " badpixlocation=PN_CCD"+pnCCD+".evts detmaptype=psf;"+ \
                           "specgroup spectrumset=PN_spectrum_source_"+obsID+".fits"+ \
                               " mincounts=20 oversample=3 rmfset=PN_"+obsID+".rmf"+ \
                               " arfset=PN_"+obsID+".arf backgndset=PN_spectrum_background_"+obsID+".fits"+ \
                               " groupedset=PN_spectrum_grouped_"+obsID+".fits;"
                           #"fv PN_spectrum_grouped_"+obsID+".fits;"
                           , shell=True)
            print('\nPN Spectra extracted.')
            print('\nExtracting Image Mode Spectra for obsID {} finished.'.format(obsID)) 

        return print('\nCompleted extracting Spectra in Image Mode.')


    ##-- fit Spectra to extracted group Spectra data --##
    def xspec_fitSpectra (self, instName='MOS1_CCD1'):
        """
        Fits a model Spectra to the extracted Spectra data using xspec package.
        See: https://www.cosmos.esa.int/web/xmm-newton/sas-thread-xspec

        Input: Extracted Spectra FITS file.
        Output: Spectra plots.

        THIS CANNOT'T RUN THROUGH THE PYTHON PROGRAM.
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

            models = ['zashift*powerlaw', 'zashift*bbody',
                      'tbabs*zashift*powerlaw', 'tbabs*zashift*bbody',
                      'zashift*\(powerlaw+bbody\)', 'tbabs*zashift*\(powerlaw+bbody\)']

            for model in models:
                subprocess.run("cd "+workdir+";"+ \
                               ". $HEADAS/headas-init.sh;"+ \
                               ". $SAS_DIR/setsas.sh;"+ \
                               #'''export SAS_CCF="`pwd`/ccf.cif";'''+ \
                               "python3 ~/Dropbox/Dheeraj@MIT_2020-21/pyXspec.py --obsID "+obsID+" --instName all " \
                                    +"--model "+model+" --showFig=no;"
                               , shell=True)
            print('\nSpectra fitting for obsID {} finished.'.format(obsID))
            
        return print('Fitted model Spectra to the extracted Spectra.')


    ##-- function to save the final results --##
    def save_results (self):
        """
        Does one tasks
                Copies the Spectra images from each obsID directory to a results 
                folder in the main directory.
        """
        print('\nLastly saving results for each obsID to a common results folder.')

        maindir = self.workdir
        #objName+".dat;"
        #-- make the results directory --#
        if not os.path.isdir(maindir+'/results'):
            subprocess.run("cd "+maindir+";"+ \
                           "mkdir results/", shell=True)
        resultdir = maindir+'/results'

        #-- iterating for all the obsIDs --#      
        for obsID in self.obsIDs:
            print('\nSaving results for obsID {}.'.format(obsID))
            workdir = maindir+'/'+obsID+'/work'
            
            #-- copy Event lists and other results --#
            spectraImgFiles = glob.glob(workdir+'/*spectra*.png')

            for file in spectraImgFiles:
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
    obj = spectra(ra=args.ra, dec=args.dec, workdir=args.workdir, \
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
    obj.extractSpectra_imageMode()
    #obj.xspec_fitSpectra()
    #obj.save_results()

    print('\nHurray! The Spectra method ran successfully.')
    if len(obj.smallMode) != 0:
        print(f'\nThe following obsIDs have Small-mode MOS data.\n{obj.smallMode}')
    if len(obj.badObs) != 0:
        print('These obsIDs were excluded from analysis: ', obj.badObs)

    return True


if __name__=="__main__":
    
    description = 'Program to reduce RGS and OM data and obtain the Spectra of XMM objects.'
    
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

