
######################################
###      The `epicObj` class       ###
######################################


import os
import sys
import subprocess
import requests
import wget

import glob
import pickle
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from datetime import datetime
from astropy.table import Table

from findOverlap import findOverlap
from plotAnal import plotAnal



##-- function to convert yes/no to bool --##
def strToBool (s, inverse=False):
    if type(s)==bool:
        if inverse:
            if s:
                return 'yes'
            else:
                return 'no'
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


##-- the XMM Object Class --##
class epicObj:

    ##-- initialize some common parameters --##
    def __init__ (self, ra, dec, workdir, sas_dir, headas, sas_ccfpath, \
                  srcCircRadius=None, bkgCircRadius=None, dSrcThreshold=70, bkgCircGap=2.5, \
                  edetectmode='chain', esp_nsplinenodes=14, gti_indiThreshold=500, gti_combThreshold=None, \
                  lcBinSize=25, binBkglc='no', saveFig=True, showFig=True, ignorePileup=False):
    
        self.ra = ra
        self.dec = dec
        self.workdir = workdir
        
        self.sas_dir = sas_dir
        self.headas = headas
        self.sas_ccfpath = sas_ccfpath
        
        self.obsIDs = list()
        self.badObs = list()      #--list of obsIDs, excluded from analysis, having total GTI below gti_combThreshold.
        self.sourceCCDs = {}      #--CCD numbers for each obsID.
        self.sourceLoc = {}       #--source location parameters for each obsID.
        self.backgroundLoc = {}   #--background circle coordinates and radius for each obsID.

        self.smallMode = {}       #--is obsID in smallMode, dict for all obsIDs.
        self.otherSources = {}    #--coords of other Sources for each obsID.
        
        self.srcCircRadius = srcCircRadius    #--radius of the source circle, arcsec changed to pixels later.
        self.bkgCircRadius = bkgCircRadius    #--radii of the background circles, arcsec.
        self.dSrcThreshold = dSrcThreshold    #--threshold distance from the source, arcsec.
        self.bkgCircGap = bkgCircGap          #--gap to have around the circles, in pixels.

        self.edetectmode = edetectmode               #--run edetect_chain or individual tasks?
        self.esp_nsplinenodes = esp_nsplinenodes     #--for esplinemap within edetect_chain, controls no. of sources detected.
        self.gti_indiThreshold = gti_indiThreshold   #--include individuals GTIs greater than this value, seconds.
        self.gti_combThreshold = gti_combThreshold   #--remove obsIDs having total GTI below this value from further analysis, seconds.
        self.lcBinSize = lcBinSize            #--binning size for the light curve, seconds.
        self.binBkglc = binBkglc              #--whether to bin background light curve or not?

        self.saveFig = strToBool(saveFig)
        self.showFig = strToBool(showFig)
        self.ignorePileup = strToBool(ignorePileup)

            #-- used for obtaining indi inst bkgCircs during Spectral analysis --#
        self.otherSourcesIndi = {}    #--other Sources coords for each obsIDs for each instrument.
        self.backgroundLocIndi = {}   #--background circle coords and radius for each obsID for each instrument.


    def _check_dir(self, path):
        """
        Checks if the given path exists and is writable. If not, it creates the directory and checks permissions.
        """
        if not os.path.isdir(path):
            os.makedirs(path, exist_ok=True)
        if not (os.access(path, os.W_OK) and os.access(path, os.X_OK)):
            print(f"Error: You do not have write/execute permissions in {path}.")
            sys.exit(1)

    
    ##-- function to find the obsIDs --##
    def findObsIDs (self):
        """
        :Input: Coordinates of the object and path to work directory.
        :Output: List of obsIDs.
        """
        ra, dec, workdir = str(self.ra), str(self.dec), self.workdir
        print('\nLooking for obsIDs at RA={} and DEC={}\nWORKDIR is set at {}'.format(ra,dec,workdir))
        
        #-- check if workdir exists --#
        self._check_dir(workdir)
            
        #-- check if browse_extract_wget.pl file exists --#
        if not os.path.isfile(workdir+"/"+"browse_extract_wget.pl"):
            print('\nbrowse_extract_wget.pl not found.')
            subprocess.run(f"cd '{workdir}';"+ \
                           "wget -q https://heasarc.gsfc.nasa.gov/FTP/heasarc/software/web_batch/browse_extract_wget.pl", shell=True)
            
            print('browse_extract_wget.pl downloaded.')
            print('\nPlease check the PERL path in the file. If required, correct the path given in first line and save the file.')
            subprocess.run(f"cd '{workdir}';"+ \
                           "gedit browse_extract_wget.pl &", shell=True)
            
        #-- download and save parts of xmmmaster table --##
        subprocess.run(f"cd '{workdir}';"+ \
                       "chmod +x browse_extract_wget.pl;"+ \
                       "./browse_extract_wget.pl table=xmmmaster position='"+ra+","+dec+ \
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
        subprocess.run(f"cd '{workdir}';"+ \
                       "cp obsIDs_list.dat obsIDs_list_"+objName+".dat;", shell=True) 
        
        return print('\nFound '+str(len(obsIDs))+' obsIDs for the object at given position.')
        

    ##-- function to sort the obsIDs_list.dat --##
    def sortObsIDs (self, objName=None):
        """
        To sort the obsIDs in `obsIDs_list.dat` by reading it into a pandas dataframe.
        """
        ra, dec, workdir = str(self.ra), str(self.dec), self.workdir
        print('\nReading the obsIDs to a dataframe and Sorting.')

        #-- offline usage --#
        if objName is not None:
            fname = workdir+'/obsIDs_list_'+objName+'.dat'
        else:
            fname = workdir+'/obsIDs_list.dat'
        
        #-- extract obsIDs from the file --#
        df = pd.read_csv(fname, sep='|', delim_whitespace=False, header=0)[:-1]  #--remove last line.
        cols = [s.strip() for s in df.columns.to_list()]  #--remove whitespace from column names.
        df.columns = cols
        
        	#-- to get the time strings --#
        #times = [datetime.strptime(t, '%Y-%m-%d %H:%M:%S') for t in df['time'] if type(t) is str]
        #df['time'] = pd.to_datetime(df['time'], format='%Y-%m-%d %H:%M:%S')

        obsIDs, times = [], []
        for i, s in enumerate(df['_Search_Offset'].fillna(0)):
            if s!=0 and float(s.split()[0].strip())<1.0:           #--remove those which are too far off.
                obsIDs.append( '0' + str(int(df['obsid'][i])) )    #--the database has a 0 at the start.
                times.append( datetime.strptime(df['time'][i], '%Y-%m-%d %H:%M:%S') )
                
        #-- create a new sorted df --#
        df1 = pd.DataFrame({'obsIDs':obsIDs, 'time':times})
        df1.sort_values(by=['time'], inplace=True)
        
        #print(df1)
        print('\nSorted the obsIDs.')
        return df1
        
    
    ##-- function to download the data --##
    def downloadData (self):
        """
        :Input: List of obsIDs for which data has to be downloaded.
        :Output: None.
        
        First it is checked is the the downloading has been done before.
        If yes, then if the size of the obsID folder is large, it's downloading is skipped.
        """        
        print('\nStarting Data Download.')
        workdir = self.workdir
        wgetRef = "wget -q -nH --no-check-certificate --cut-dirs=4 -r -l0 -c -N -np -R 'index*'  -erobots=off --retr-symlinks "
        wgetRef = wgetRef + "--show-progress --progress=bar:force "  #--to show progress bar.
        
        #-- iterate for all the obsIDs --#
        for obsID in self.obsIDs:
        
            downPath = workdir+'/'+obsID+'/'
                #-- check if already downloaded --#
            if os.path.isdir(downPath)==True:
                alreadyDownSize = sum(d.stat().st_size for d in os.scandir(downPath) if d.is_dir())
                print('\nFolder for obsID {} already exists with file size {}.'.format(obsID, alreadyDownSize))
                if alreadyDownSize>10000:
                    print('\n\tSkipping download for obsID {}.\n'.format(obsID))
                    continue
            
            #-- download the data --#
            url = "https://heasarc.gsfc.nasa.gov/FTP/xmm/data/rev0//"+obsID+"/."
            wgetFull = wgetRef + url
            print('\nDownloading data for obsID {} using wget.'.format(obsID))
            subprocess.run(f"cd '{workdir}';"+wgetFull, shell=True)
            
            #-- unzip downloaded files --#
            print('\nUnzipping the downloaded ODF tar files.')
            odfPath = downPath+'ODF'
            subprocess.run(f"cd '{odfPath}';"+
                           "gunzip *.gz", shell=True)
            print('Data download for obsID {} finished.'.format(obsID))
                
        return print('\nCompleted Downloading Data!\n')
        
        
    ##-- function to reduce EPIC data --##
    def reduceEPICdata (self):
        """
        Runs the XMMSAS commands on the shell to reduce EPIC data.
        
        :Input: SAS Summary file in the work directory, that is downloaded alongwith the data.
        :Output: EPIC PN and MOS1&2 calibrated and concatenated Event lists
                produced by SAS tasks `epproc` and `emproc`.
                
        These event lists contain Instrumental GTI for each of the CCDs.
        See Notes and http://xmm-tools.cosmos.esa.int/external/sas/current/doc/epicproc/node17.html
        
        This function largely follows the following SAS thread:
        https://www.cosmos.esa.int/web/xmm-newton/sas-thread-epic-reprocessing
        """
        print('\nInitiating XMMSAS commands for EPIC Data Reduction.')
        workdir = self.workdir
        
        #-- set the environment variables --#
        os.environ['SAS_DIR'] = self.sas_dir
        os.environ['HEADAS'] = self.headas
        os.environ['SAS_CCFPATH'] = self.sas_ccfpath
        
        #-- run the XMMSAS commands --#
        for obsID in self.obsIDs:
            print('\nReducing data for obsID {}.'.format(obsID))
            os.environ['SAS_ODF'] = workdir+'/'+obsID+'/ODF'
            
            savedir = workdir+'/'+obsID+'/work'
            self._check_dir(savedir)
                
            subprocess.run(f"cd '{savedir}';"+ \
                           ". $HEADAS/headas-init.sh;"+ \
                           ". $SAS_DIR/setsas.sh;"+ \
                           #sasversion;+ \
                           "cifbuild;"+ \
                           '''export SAS_CCF="`pwd`/ccf.cif";'''+ \
                           #"echo $SAS_CCF;"+ \
                           "odfingest;"+ \
                           '''export SAS_ODF="`pwd`/` ls *SUM.SAS`";'''+ \
                           #"echo $SAS_ODF;"+ \
                           "epproc;"+ \
                           "emproc;", shell=True)
            print('Data reduction for obsID {} finished.'.format(obsID))
        
        return print('\nCompleted reducing the data and created EPIC Event lists!')


    ##-- function to extract Flare GTI --##
    def extract_flareGTI (self):
        """
        Runs the XMMSAS commands to extract Flare GTI from PN Event List.
        
        :Input: Unfiltered EPIC PN and MOS1&2 Events lists.
        :Output: `flareGTI` file for the obsIDs.
        
        This function largely follows the following SAS thread:
        https://www.cosmos.esa.int/web/xmm-newton/sas-thread-epic-filterbackground
        
        Only one of either EPIC PN or MOS Flare Background filtering is required 
        since the Flare affects the whole telescope.
        """
        print('\nStarting EPIC data filtering to extract Flare GTIs.')
        
        #-- set the environment variables --#
        os.environ['SAS_DIR'] = self.sas_dir
        os.environ['HEADAS'] = self.headas
        os.environ['SAS_CCFPATH'] = self.sas_ccfpath
        
        #-- iterating for all the obsIDs --#
        for obsID in self.obsIDs:
            print('\nRunning commands for obsID {}.'.format(obsID))
            workdir = self.workdir+'/'+obsID+'/work'
            
            #-- grab file names --#
            pnFile = glob.glob(workdir+'/*EPN*ImagingEvts*')[0]
            
            #-- filter EPIC PN data --#
            subprocess.run(f"cd '{workdir}';"+ \
                           ". $HEADAS/headas-init.sh;"+ \
                           ". $SAS_DIR/setsas.sh;"+ \
                           #"fv "+pnFile+";"+ \
                           #"sasversion;"+ \
                           "evselect table="+pnFile+ \
                               " withrateset=Y rateset=ratePN.fits maketimecolumn=Y timebinsize=100 \
                               makeratecolumn=Y expression='#XMMEA_EP && (PI>10000&&PI<12000) && (PATTERN==0)';"+ \
                           #"dsplot table=ratePN.fits x=TIME y=RATE.ERROR;"+ \
                           "tabgtigen table=ratePN.fits expression='RATE<0.4' gtiset=flareGTI.fits;"
                           #"fv flareGTI.fits;"
                           , shell=True)
            print('\nFlare GTI extracted for obsID {} finished.'.format(obsID))

        return print('\nExtracted the Flare GTIs!')
        

    ##-- function to remove Flare background --##
    def removeFlareBackground (self):
        """
        Runs the XMMSAS commands to remove Flare background from EPIC Event List
        and obtained cleaned event lists for each EPIC instrument.
        The Flare GTI file extracted previously using `extract_flareGTI` is used.
        Since the Flare Background is same for all the instruments onboard, 
        this same `flareGTI` file can be used for them all.

        Note that these Flare Background filtered Event Lists are only used for 
        Spectral Analysis, since for Light Curve extraction, `combinedGTI` file 
        and combined `PNMOS12.evts` EPIC data is used.
        
        :Input: Unfiltered EPIC PN and MOS1&2 Events lists and `flareGTI` file.
        :Output: Flare Background filtered Event Lists for each EPIC instrument.
        
        This function largely follows the following SAS thread:
        https://www.cosmos.esa.int/web/xmm-newton/sas-thread-epic-filterbackground
        """
        print('\nStarting EPIC data filtering to get Flare GTIs.')
        
        #-- set the environment variables --#
        os.environ['SAS_DIR'] = self.sas_dir
        os.environ['HEADAS'] = self.headas
        os.environ['SAS_CCFPATH'] = self.sas_ccfpath
        
        #-- iterating for all the obsIDs --#
        for obsID in self.obsIDs:
            print('\nRunning commands for obsID {}.'.format(obsID))
            workdir = self.workdir+'/'+obsID+'/work'
            
            #-- grab file names --#
            pnFile = glob.glob(workdir+'/*EPN*ImagingEvts*')[0]
            mos1File = glob.glob(workdir+'/*EMOS1*ImagingEvts*')[0]
            mos2File = glob.glob(workdir+'/*EMOS2*ImagingEvts*')[0]
            
            #-- filter EPIC PN data --#
            subprocess.run(f"cd '{workdir}';"+ \
                           ". $HEADAS/headas-init.sh;"+ \
                           ". $SAS_DIR/setsas.sh;"+ \
                           "evselect table="+pnFile+ \
                               " withfilteredset=Y filteredset=PNclean.ds destruct=Y keepfilteroutput=T \
                               expression='#XMMEA_EP && gti(flareGTI.fits, TIME) && (PI>150)';"+ \
                           #"fv PNclean.ds;"+ \
                           "evselect table="+mos1File+ \
                               " withfilteredset=Y filteredset=MOS1clean.ds destruct=Y keepfilteroutput=T \
                               expression='#XMMEA_EM && gti(flareGTI.fits, TIME) && (PI>150)';"+ \
                           #"fv PNclean.ds;"+ \
                           "evselect table="+mos2File+ \
                               " withfilteredset=Y filteredset=MOS2clean.ds destruct=Y keepfilteroutput=T \
                               expression='#XMMEA_EM && gti(flareGTI.fits, TIME) && (PI>150)';"
                           #"fv PNclean.ds;"+ \
                           , shell=True)
            print('\nFlare GTI extracted for obsID {} finished.'.format(obsID))

        return print('\nCompleted EPIC PN data filtering and extracted Flare GTIs!')
        
        
    ##-- function to find where Source is located --##
    def findSourceCCD (self):
        """
        Uses 'ecoordconv' command of SAS to find which CCD contains the Source
        in each of the EPIC's three cameras.
        
        :Input: EPIC concatenated and calibrated Event lists obtained by running 
               emproc and epproc.
        :Output: None.
        
        self.sourceCCDs is updated as a dictionary of dictionaries containing the
        Source CCD number for PN, MOS1 and MOS2 cameras for all the obsIDs.
        """
        print('\nFinding which CCDs contain the Source.')
        
        #-- set the environment variables --#
        os.environ['SAS_DIR'] = self.sas_dir
        os.environ['HEADAS'] = self.headas
        os.environ['SAS_CCFPATH'] = self.sas_ccfpath
        
        ra, dec = str(self.ra), str(self.dec)
        
        #-- iterating for all the obsIDs --#
        for obsID in self.obsIDs:
            print('\nFinding Source CCD for obsID {}.'.format(obsID))
            workdir = self.workdir+'/'+obsID+'/work'
            ccdList = {'PN':None, 'MOS1':None, 'MOS2':None}
            locParams = {'x':None, 'y':None, 'r':None}
            
            #-- grab file names --#
            pnFile = glob.glob(workdir+'/*EPN*ImagingEvts*')[0]
            mos1File = glob.glob(workdir+'/*EMOS1*ImagingEvts*')[0]
            mos2File = glob.glob(workdir+'/*EMOS2*ImagingEvts*')[0]
            
            #-- find Source CCD in EPIC PN data --#
            cmdOut = subprocess.run(f"cd '{workdir}';"+ \
                           ". $HEADAS/headas-init.sh;"+ \
                           ". $SAS_DIR/setsas.sh;"+ \
                           '''export SAS_CCF="`pwd`/ccf.cif";'''+ \
                           "ecoordconv imageset="+pnFile+":EVENTS srcexp='' \
                               x="+ra+" y="+dec+" pos2eqpos=no im2eqpos=no coordtype=eqpos;", \
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            line = cmdOut.stdout.splitlines()[-3]
            pnCCD = int(line.split()[-1])
            ccdList['PN'] = pnCCD
            print('\nSource CCD in PN:', pnCCD)
            
            #-- find Source CCD in EPIC MOS1 data --#
            cmdOut = subprocess.run(f"cd '{workdir}';"+ \
                           ". $HEADAS/headas-init.sh;"+ \
                           ". $SAS_DIR/setsas.sh;"+ \
                           '''export SAS_CCF="`pwd`/ccf.cif";'''+ \
                           "ecoordconv imageset="+mos1File+":EVENTS srcexp='' \
                               x="+ra+" y="+dec+" pos2eqpos=no im2eqpos=no coordtype=eqpos;", \
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            line = cmdOut.stdout.splitlines()[-3]
            mos1CCD = int(line.split()[-1])
            ccdList['MOS1'] = mos1CCD
            print('Source CCD in MOS1:', mos1CCD)
            
            #-- find Source CCD in EPIC MOS2 data --#
            cmdOut = subprocess.run(f"cd '{workdir}';"+ \
                           ". $HEADAS/headas-init.sh;"+ \
                           ". $SAS_DIR/setsas.sh;"+ \
                           '''export SAS_CCF="`pwd`/ccf.cif";'''+ \
                           "ecoordconv imageset="+mos2File+":EVENTS srcexp='' \
                               x="+ra+" y="+dec+" pos2eqpos=no im2eqpos=no coordtype=eqpos;", \
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            line = cmdOut.stdout.splitlines()[-3]
            mos2CCD = int(line.split()[-1])
            ccdList['MOS2'] = mos2CCD
            print('Source CCD in MOS2:', mos2CCD)
            
            #-- find Source X and Y --#
            line2 = cmdOut.stdout.splitlines()[-8]
            x, y = float(line2.split()[2]), float(line2.split()[3])
            locParams['x'], locParams['y'] = x, y
            print('Source X and Y in Sky coords: {}, {}'.format(x, y))
            
            #-- save ccdList and sourceLoc --#
            self.sourceCCDs[obsID] = ccdList
            self.sourceLoc[obsID] = locParams
            print('\nFinding Source CCD for obsID {} finished.'.format(obsID))
        
        #print(self.sourceLoc)
        return print('\nFound the CCDs containing the Source!')
        
    
    ##-- function to extract Instrumental GTIs --##
    def extract_InstrumentalGTIs (self):
        """
        Extracts Instrumental GTIs from EPIC data using 'fextract' command.
        
        :Input: EPIC concatenated and calibrated Event lists obtained by running 
               emproc and epproc.
        :Output: Intrumental GTI FITS files for EPIC PN, MOS1 and MOS2.
        
        'fextract' command can't overwrite the extracted file if it already exists 
        in the work directory, so they are deleted using 'rm -rf', which doesn't 
        raises any error even if no such files are found for deleting during the 
        first run.
        """
        print('\nStarting Instrumental GTIs Extraction.')
        
        #-- set the environment variables --#
        os.environ['SAS_DIR'] = self.sas_dir
        os.environ['HEADAS'] = self.headas
        os.environ['SAS_CCFPATH'] = self.sas_ccfpath
        
        #-- iterating for all the obsIDs --#      
        for obsID in self.obsIDs:
            print('\nGetting Instrumental GTIs for obsID {}.'.format(obsID))
            workdir = self.workdir+'/'+obsID+'/work'
            
            #-- remove previously extracted files --#
            subprocess.run(f"cd '{workdir}';"+ \
                           "rm -rf pnGTI.fits;"+ \
                           "rm -rf mos1GTI.fits;"+ \
                           "rm -rf mos2GTI.fits;", shell=True)
            
            #-- grab file names --#
            pnFile = glob.glob(workdir+'/*EPN*ImagingEvts*')[0]
            mos1File = glob.glob(workdir+'/*EMOS1*ImagingEvts*')[0]
            mos2File = glob.glob(workdir+'/*EMOS2*ImagingEvts*')[0]
            
            #-- Source CCD numbers --#
            pnCCD = str(self.sourceCCDs[obsID]['PN'])
            mos1CCD = str(self.sourceCCDs[obsID]['MOS1'])
            mos2CCD = str(self.sourceCCDs[obsID]['MOS2'])
            
            #-- extract EPIC PN, MOS1 and MOS2 data --#
            subprocess.run(f"cd '{workdir}';"+ \
                           ". $HEADAS/headas-init.sh;"+ \
                           ". $SAS_DIR/setsas.sh;"+ \
                           #"fv "+pnFile+";"+ \
                           "fextract "+pnFile+"[STDGTI"+pnCCD.zfill(2)+"] "+workdir+"/pnGTI.fits;"
                           #"fv pnGTI.fits"+ \
                           #"fv "+mos1File+";"+ \
                           "fextract "+mos1File+"[STDGTI"+mos1CCD.zfill(2)+"] "+workdir+"/mos1GTI.fits;"
                           #"fv mos1GTI.fits"+ \
                           #"fv "+mos2File+";"+ \
                           "fextract "+mos2File+"[STDGTI"+mos2CCD.zfill(2)+"] "+workdir+"/mos2GTI.fits;"
                           #"fv mos2GTI.fits;"
                           , shell=True)
            print('Instrumental GTIs for obsID {} obtained for EPIC PN, MOS1 and MOS2 data.'.format(obsID))
        
        return print('\nSuccessfully extracted Instrumental GTIs!')


    ##-- function to combine GTIs --##
    def combineGTIs (self):
        """
        Combines Flare and Instrumental GTIs.
        
        :Input: Flare, EPIC PN, MOS1 and MOS2 GTIs.
        :Output: Combined Flare and Instrumental GTI FITS file.
        """
        print('\nStarting process to Combine Flare, EPIC PN, MOS1 and MOS2 GTIs.')
        
        #-- set the environment variables --#
        os.environ['SAS_DIR'] = self.sas_dir
        os.environ['HEADAS'] = self.headas
        os.environ['SAS_CCFPATH'] = self.sas_ccfpath

        #-- combined GTI threshold --#
        if self.gti_combThreshold==None:
            gti_combThreshold = 2*float(self.gti_indiThreshold)
        else:
            gti_combThreshold = float(self.gti_combThreshold)
        print(f'gti_combThreshold is set at {gti_combThreshold}')

        goodObs = []  #--list of obsIDs with GTI total more than gti_combThreshold.
        
        #-- iterating for all the obsIDs --#      
        for obsID in self.obsIDs:
            print('\nCombining GTIs for obsID {}.'.format(obsID))
            workdir = self.workdir+'/'+obsID+'/work'
            
            #-- remove previously combined data files --#
            subprocess.run(f"cd '{workdir}';"+ \
                           "rm -rf combinedGTI_"+obsID+".fits;", shell=True)
                           
            #-- combine the GTIs --#
            subprocess.run(f"cd '{workdir}';"+ \
                           ". $HEADAS/headas-init.sh;"+ \
                           ". $SAS_DIR/setsas.sh;"+ \
                           "mgtime "+workdir+"/pnGTI.fits,"+workdir+"/mos1GTI.fits,"+workdir+"/mos2GTI.fits,"+ \
                               workdir+"/flareGTI.fits "+workdir+"/combinedGTI_"+obsID+".fits AND;"
                           #"ftmgtime "+workdir+"/pnGTI.fits,"+workdir+"/mos1GTI.fits"+workdir+"/mos2GTI.fits "+ \
                           #    workdir+"/flareGTI.fits "+workdir+"/combinedGTI_"+obsID+".fits AND;"
                           #"fv combinedGTI_"+obsID+".fits;"
                           , shell=True)

            #-- check if combined GTI length is less --#
            print('Checking length of combined GTI.')
            tgti = Table.read(workdir+'/'+'combinedGTI_'+obsID+'.fits', format='fits', hdu=1)
            tgti.write(workdir+'/'+'combinedGTI_'+obsID+'.csv', overwrite=True)
            dfgti = pd.read_csv(workdir+'/'+'combinedGTI_'+obsID+'.csv')
            gtiStart, gtiStop = dfgti.START, dfgti.STOP

            lengths = [stop-start for start, stop in zip(gtiStart, gtiStop)]
            gtiTotal = sum(lengths)
            if gtiTotal >= gti_combThreshold:
                goodObs.append(obsID)
                print(f'Total GTI length is {gtiTotal} seconds.')
            else:
                print(f'Total GTI length for obsID {obsID} is below threshold. ObsID would be excluded from further analysis.')

            print('\nCombining GTIs for obsID {} finished.'.format(obsID))

        #-- save good and bad obsIDs --#
        badObs = list( set(self.obsIDs) - set(goodObs) )
        self.obsIDs = goodObs
        self.badObs = badObs
            
        return print('\nSuccessfully combined Flare and Instrumental GTIs!')
       
        
    ##-- function to combine EPIC data --##
    def combineEPICdata (self):
        """
        Combines EPIC PN, MOS1 and MOS2 Event lists data corresponding to the CCD
        where the Source is located.
        
        :Input: EPIC concatenated and calibrated Event lists obtained by running 
               emproc and epproc.
        :Output: Combined PNMOS12.evts and MOS12.evts Event lists.
        """
        print('\nStarting process to Combine EPIC PN, MOS1 and MOS2 data.')
        
        #-- set the environment variables --#
        os.environ['SAS_DIR'] = self.sas_dir
        os.environ['HEADAS'] = self.headas
        os.environ['SAS_CCFPATH'] = self.sas_ccfpath
        
        #-- iterating for all the obsIDs --#      
        for obsID in self.obsIDs:
            print('\nCombining EPIC data for obsID {}.'.format(obsID))
            workdir = self.workdir+'/'+obsID+'/work'
            
            #-- remove previously combined data files --#
            subprocess.run(f"cd '{workdir}';"+ \
                           "rm -rf MOS12.evts;"+ \
                           "rm -rf PNMOS12.evts;", shell=True)
            
            #-- grab file names --#
            pnFile = glob.glob(workdir+'/*EPN*ImagingEvts*')[0]
            mos1File = glob.glob(workdir+'/*EMOS1*ImagingEvts*')[0]
            mos2File = glob.glob(workdir+'/*EMOS2*ImagingEvts*')[0]
            
            #-- Source CCD numbers --#
            pnCCD = str(self.sourceCCDs[obsID]['PN'])
            mos1CCD = str(self.sourceCCDs[obsID]['MOS1'])
            mos2CCD = str(self.sourceCCDs[obsID]['MOS2'])
            
            #-- extract corresponding CCD data --#
            subprocess.run(f"cd '{workdir}';"+ \
                           ". $HEADAS/headas-init.sh;"+ \
                           ". $SAS_DIR/setsas.sh;"+ \
                           "evselect table="+pnFile+ \
                               " withfilteredset=yes filteredset="+workdir+"/PN_CCD"+pnCCD+".evts"+ \
                               " keepfilteroutput=yes expression='#XMMEA_EP && (CCDNR=="+pnCCD+") && (PATTERN<=4)';"+ \
                           #"fv PN_CCD"+pnCCD+".evts;"+ \
                           "evselect table="+mos1File+ \
                               " withfilteredset=yes filteredset="+workdir+"/MOS1_CCD"+mos1CCD+".evts"+ \
                               " keepfilteroutput=yes expression='#XMMEA_EM && (CCDNR=="+mos1CCD+") && (PATTERN<=12)';"+ \
                           "evselect table="+mos2File+ \
                               " withfilteredset=yes filteredset="+workdir+"/MOS2_CCD"+mos2CCD+".evts"+ \
                               " keepfilteroutput=yes expression='#XMMEA_EM && (CCDNR=="+mos2CCD+") && (PATTERN<=12)';", shell=True)
            print('\nEPIC PN, MOS1 and MOS2 data from corresponding CCD extracted.')
            
            #-- combine MOS1 and MOS2 data --#
            subprocess.run(f"cd '{workdir}';"+ \
                           ". $HEADAS/headas-init.sh;"+ \
                           ". $SAS_DIR/setsas.sh;"+ \
                           "merge set1="+workdir+"/MOS1_CCD"+mos1CCD+".evts set2="+workdir+"/MOS2_CCD"+mos2CCD+".evts"+ \
                           " outset="+workdir+"/MOS12.evts", shell=True)
            print('\nMOS1_CCD and MOS2_CCD data combined.')
                        
            #-- combine MOS12 and PN data --#
            subprocess.run(f"cd '{workdir}';"+ \
                           ". $HEADAS/headas-init.sh;"+ \
                           ". $SAS_DIR/setsas.sh;"+ \
                           "merge set1="+workdir+"/MOS12.evts set2="+workdir+"/PN_CCD"+pnCCD+".evts"+ \
                           " outset="+workdir+"/PNMOS12.evts", shell=True)
            print('\nMOS12 and PN_CCD data combined.')
                           
            print('\nCombining data for obsID {} finished.'.format(obsID))
            
        return print('\nSuccessfully combined EPIC PN and MOS data!')


    ##-- to find other sources --##
    def find_otherSources (self):
        """
        Function to find all the other Sources, in addition to the main Source, lying within the region of interest.
        See: https://www.cosmos.esa.int/web/xmm-newton/sas-thread-src-find-stepbystep
    
        :Input: PNMOS12.evts and PN_CCD##.evts and corresponding extracted images in the Energy band 0.3-10 KeV
               using evselect.
        :Output: Exposer maps, masks, the csv file containing the x and y coordinates of the other Sources, 
                apart from the region file from DS9.

        Note that other sources may be found within any of the regions using different input files,
        for instance 'PN_image_full.fits' obtained from 'PNclean.ds' etc.

        However, only the sources within the PNMOS12 overlap region or the PN region, in case of small mode,
        are of interest. Since, the PN region necessarily contains the overlap region within it regardless of 
        the mode of observation, this should have been sufficient to get other sources. However, the MOS data adds 
        some of its own sources and 'edetect_chain' is easier to use and can be used with multiple imagesets,
        so both 'PNMOS12' and 'PN_CCD##' are being used to get the sources.

        Further, right now the SAS command 'edetect_chain', which is pretty slow when number of sources is large,
        is been used. Individual tasks within 'edetect_chain' can be run separately to fasten this up, although
        that would require changes in the createOutputImages() function below.
        """
        print('\nFinding all the other Sources.')
        
        #-- set the environment variables --#
        os.environ['SAS_DIR'] = self.sas_dir
        os.environ['HEADAS'] = self.headas
        os.environ['SAS_CCFPATH'] = self.sas_ccfpath
    
        esp_nsplinenodes = str(self.esp_nsplinenodes)

        #-- iterating for all the obsIDs --#     
        for obsID in self.obsIDs:
            print('\nFinding other Sources for obsID {}.'.format(obsID))
            workdir = self.workdir+'/'+obsID+'/work'

            pnCCD = str(self.sourceCCDs[obsID]['PN'])

            #-- grab Attitude File --#
            AttFile = glob.glob(workdir+'/*AttHk*.ds')[0]
            """
            subprocess.run(f"cd '{workdir}';"+ \
                           ". $HEADAS/headas-init.sh;"+ \
                           ". $SAS_DIR/setsas.sh;"+ \
                           '''export SAS_CCF="`pwd`/ccf.cif";'''+ \
                           "evselect table=PN_CCD"+pnCCD+".evts:EVENTS \
                               imagebinning='binSize' imageset='PN_CCD"+pnCCD+"_image.fits' withimageset=yes \
                               xcolumn='X' ycolumn='Y' ximagebinsize=80 yimagebinsize=80 \
                               expression='#XMMEA_EP && (PI in [300:10000]) && (PATTERN in [0:12])';"+ \
                           #"ds9 PN_CCD"+pnCCD+"_image.fits -scale log -cmap bb;"
                           "eexpmap attitudeset="+AttFile+" eventset=PN_CCD"+pnCCD+".evts imageset=PN_CCD"+pnCCD+"_image.fits \
                               expimageset=PN_CCD"+pnCCD+"_expmap.ds pimin='300' pimax='10000';"+ \
                           #"fv PN_CCD"+pnCCD+"_expmap.ds;"
                           "emask expimageset=PN_CCD"+pnCCD+"_expmap.ds threshold1=0.25 detmaskset=PN_CCD"+pnCCD+"_mask.ds;"+ \
                           #"fv PN_CCD"+pnCCD+"_mask.ds;"
                           "eboxdetect usemap=no likemin=8 withdetmask=yes detmasksets=PN_CCD"+pnCCD+"_mask.ds \
                               imagesets=PN_CCD"+pnCCD+"_image.fits expimagesets=PN_CCD"+pnCCD+"_expmap.ds \
                               pimin=300 pimax=10000 boxlistset=eboxlist_local.fits;"+ \
                           "esplinemap bkgimageset=PN_CCD"+pnCCD+"_bkg.ds scut=0.005 imageset=PN_CCD"+pnCCD+"_image.fits \
                               nsplinenodes=16 withdetmask=yes detmaskset=PN_CCD"+pnCCD+"_mask.ds withexpimage=yes \
                               expimageset=PN_CCD"+pnCCD+"_expmap.ds boxlistset=eboxlist_local.fits;" + \
                           #"fv PN_CCD"+pnCCD+"_bkg.ds;" 
                           "eboxdetect usemap=yes bkgimagesets=PN_CCD"+pnCCD+"_bkg.ds likemin=8 withdetmask=yes \
                               detmasksets=PN_CCD"+pnCCD+"_mask.ds imagesets=PN_CCD"+pnCCD+"_image.fits \
                               expimagesets=PN_CCD"+pnCCD+"_expmap.ds pimin=300 pimax=10000 boxlistset=eboxlist_map.fits;"+ \
                           "emldetect imagesets=PN_CCD"+pnCCD+"_image.fits expimagesets=PN_CCD"+pnCCD+"_expmap.ds \
                               bkgimagesets=PN_CCD"+pnCCD+"_bkg.ds boxlistset=eboxlist_map.fits ecf=2.0 \
                               mllistset=emllist.fits mlmin=10 determineerrors=yes;" + \
                           "esensmap expimagesets=PN_CCD"+pnCCD+"_expmap.ds bkgimagesets=PN_CCD"+pnCCD+"_bkg.ds \
                               detmasksets=PN_CCD"+pnCCD+"_mask.ds mlmin=10 sensimageset=PN_CCD"+pnCCD+"_sens_map.fits;"+ \
                           #"ds9 PN_CCD"+pnCCD+"_sens_map.fits -scale log -cmap bb;" 
                           "srcdisplay boxlistset=emllist.fits imageset=PN_CCD"+pnCCD+"_image.fits sourceradius=0.005 \
                               withregionfile=true regionfile=allSources.reg;"
                           , shell=True)
            """
            subprocess.run(f"cd '{workdir}';"+ \
                           ". $HEADAS/headas-init.sh;"+ \
                           ". $SAS_DIR/setsas.sh;"+ \
                           '''export SAS_CCF="`pwd`/ccf.cif";'''+ \
                           "evselect table=PNMOS12.evts:EVENTS imagebinning='binSize' \
                               imageset='PNMOS12_image_full.fits' withimageset=yes \
                               xcolumn='X' ycolumn='Y' ximagebinsize=80 yimagebinsize=80 \
                               expression='(PI in [300:10000]) && (PATTERN in [0:12])';"+ \
                           "evselect table=PN_CCD"+pnCCD+".evts:EVENTS \
                               imagebinning='binSize' imageset='PN_CCD"+pnCCD+"_image.fits' withimageset=yes \
                               xcolumn='X' ycolumn='Y' ximagebinsize=80 yimagebinsize=80 \
                               expression='#XMMEA_EP && (PI in [300:10000]) && (PATTERN in [0:12])';"+ \
                           "edetect_chain imagesets='PNMOS12_image_full.fits PN_CCD"+pnCCD+"_image.fits' \
                               eventsets='PNMOS12.evts PN_CCD"+pnCCD+".evts' attitudeset="+AttFile+" \
                               pimin='300 300' pimax='10000 10000' ecf='2.0 2.0' \
                               esp_nsplinenodes="+esp_nsplinenodes+" esen_mlmin=10;"+ \
                           "srcdisplay boxlistset=emllist.fits imageset=PN_CCD"+pnCCD+"_imagesmap.fits sourceradius=0.005 \
                               withregionfile=true regionfile=allSources.reg;"
                           #"fv emllist.fits;" 
                           #"ds9 PNMOS12_image_fullsmap.fits -regions load allSources.reg \
                           #    -cmap bb -scale log -zoom 4 -export 'allSources.jpeg' 300;"
                           , shell=True) 

            t = Table.read(workdir+'/'+'emllist.fits', format='fits')
            #t = t[list(t.columns)]
            t.write(workdir+'/'+'emllist.csv', overwrite=True)
            df = pd.read_csv(workdir+'/'+'emllist.csv').drop_duplicates(['X_IMA'])
            allX, allY = np.array(df.X_IMA), np.array(df.Y_IMA)

            self.otherSources[obsID] = [(x, y) for x, y in zip(allX, allY)]

            print('\nAll other Sources for obsID {} found.'.format(obsID))
        #print(self.otherSources)    
        return print('\nFinished finding all the other Sources!')

            
    ##-- function to get Background Circles --##
    def getBackgroundCircles (self):
        """
        Calls `findOverlap.py` functions to automatically detect the overlap PNMOS12 region and 
        find 2 background cicles in it. 
        If the MOS mode is `small`, then the PN image is used for getting the background circles.

        :Input: PN_image.fits, MOS12_image.fits and PNMOS12_image.fits to get the overlap from.
        :Output: x, y coords of the background circles alongwith the radii of these circles and the
                source circle. Units of radii are pixels or Sky Coord units.
                Updates the ``sourceLoc`` and ``backgroundLoc`` dictionaries with the location
                parameters.
        """
        print('\nAutomatically detecting the Background Circles.')
        
        #-- set the environment variables --#
        os.environ['SAS_DIR'] = self.sas_dir
        os.environ['HEADAS'] = self.headas
        os.environ['SAS_CCFPATH'] = self.sas_ccfpath
        
        #-- iterating for all the obsIDs --#      
        for obsID in self.obsIDs:
            print('\nDetecting the PNMOS12 overlap region for obsID {}.'.format(obsID))
            workdir = self.workdir+'/'+obsID+'/work'
            bLocParams = {'Bx1':0, 'By1':0, 'Br1':0, 'Bx2':0, 'By2':0, 'Br2':0}
            
            #-- the source location --#
            pnCCD = str(self.sourceCCDs[obsID]['PN'])
            locParams = self.sourceLoc[obsID]
            srcX, srcY = locParams['x'], locParams['y']
            srcCoords = (srcX, srcY)

            #-- coordinates of other sources --#
            otherSrc = self.otherSources[obsID]
            
            #-- find the overlap region and background circles --#
            findOverlapObj = findOverlap(workdir=workdir, obsID=obsID, pnCCD=pnCCD, srcCoords=srcCoords, otherSrc=otherSrc, \
                                         srcR=self.srcCircRadius, bkgR=self.bkgCircRadius, srcThreshold=self.dSrcThreshold, \
                                         gap=self.bkgCircGap, saveFig=self.saveFig, showFig=self.showFig)
            bCircles, correctSrc, srcR, isSmallMode = findOverlapObj.main()
            
            #-- save the background circle parameters --#
            for i, bCircle in enumerate(bCircles):
                xKey, yKey, rKey = [key+str(i+1) for key in ['Bx', 'By', 'Br']]
                bLocParams[xKey], bLocParams[yKey], bLocParams[rKey] = bCircle[0], bCircle[1], bCircle[2]
            self.backgroundLoc[obsID] = bLocParams

            #-- save the corrected main source coordinates --#
            self.sourceLoc[obsID]['x'], self.sourceLoc[obsID]['y'] = correctSrc[0], correctSrc[1]
            self.sourceLoc[obsID]['r'] = srcR
            self.smallMode[obsID] = isSmallMode
            
            print('Background for obsID {} found.'.format(obsID))
          
        #print(self.backgroundLoc)
        return print('\nCompleted finding the Background Circles!')


    ##-- to find other sources --##
    def find_otherSources_indi (self, inst='PN'):
        """
        Function to find all the other Sources, in addition to the main Source, lying within the region of interest.
        See: https://www.cosmos.esa.int/web/xmm-newton/sas-thread-src-find-stepbystep
    
        Args:
            inst (str): name of the instrument of use.

        :Input: `inst_CCD##.evts` Event list.

        Raises:
            Error: if `inst` is not in ['PN', 'MOS1', 'MOS2']

        :Output:
            Exposer maps, masks, the csv file containing the x and y coordinates of the other Sources, 
            apart from the region file from DS9.
            Also extracts `inst_image_CCD##.fits` file which is later used to find srcBkg circs.
            Updates the ``otherSourceIndi`` dictionary for each of the obsIDs in ``self.obsIDs``

        NOTES
        -----
        This separate function based on `find_otherSources` is meant to be used for
        finding individual instrument other sources. It uses ``edetect_chain`` function 
        from XMMSAS for detecting the Sources.
        `inst_CCD##.evts` is required to extract corresponding images in the Energy band 
        0.3-10 KeV using ``evselect``.

        * The indi inst files have ``inst`` appended at the end of their names.
        * `emllist.fits` originally found for the overlap is rewritten with every run
          of ``edetect_chain``. However, the `emllist.csv` is not b'cuz of the above point
          So `emllist.csv` should be used for the overlap cases.
          And anyway, another run of the original `find_otherSources` will update the said file.

        TODO
        ----
        Check is using ``PATTERN in [0:12]`` in ``evselect expression`` for PN is correct?
        For data extraction elsewhere, the expression for PN reads, ``(PATTERN <= 4)``
        """
        print('\nFinding all the other Sources.')
        
        #-- set the environment variables --#
        os.environ['SAS_DIR'] = self.sas_dir
        os.environ['HEADAS'] = self.headas
        os.environ['SAS_CCFPATH'] = self.sas_ccfpath
    
        esp_nsplinenodes = str(self.esp_nsplinenodes)

        #-- iterating for all the obsIDs --#     
        for obsID in self.obsIDs:
            print('\nFinding other Sources for obsID {}.'.format(obsID))
            workdir = self.workdir+'/'+obsID+'/work'

            ccd = str(self.sourceCCDs[obsID][inst])
            name = inst + '_CCD' + ccd

            self.otherSourcesIndi[obsID] = {}   #--initialize empty dict for saving indi inst other Sources.

            evselectForPNMOS = "evselect table=PNMOS12.evts:EVENTS imagebinning='binSize' \
                               imageset='PNMOS12_image_full.fits' withimageset=yes \
                               xcolumn='X' ycolumn='Y' ximagebinsize=80 yimagebinsize=80 \
                               expression='(PI in [300:10000]) && (PATTERN in [0:12])';"

            pimin, pimax, ecf = '300', '10000', '2.0'

            if inst == 'PN':
                expression = "#XMMEA_EP && (PI in [300:10000]) && (PATTERN in [0:12])"
            else:
                expression = "#XMMEA_EM && (PI in [300:10000]) && (PATTERN in [0:12])"

            #-- grab Attitude File --#
            AttFile = glob.glob(workdir+'/*AttHk*.ds')[0]

            #-- run the functions --#
            subprocess.run(f"cd '{workdir}';"+ \
                           ". $HEADAS/headas-init.sh;"+ \
                           ". $SAS_DIR/setsas.sh;"+ \
                           '''export SAS_CCF="`pwd`/ccf.cif";'''+ \
                           "evselect table="+name+".evts:EVENTS \
                               imagebinning='binSize' imageset='"+name+"_image.fits' withimageset=yes \
                               xcolumn='X' ycolumn='Y' ximagebinsize=80 yimagebinsize=80 \
                               expression='"+expression+"';"+ \
                           "edetect_chain imagesets='"+name+"_image.fits' \
                               eventsets='"+name+".evts' attitudeset="+AttFile+" \
                               pimin='"+pimin+"' pimax='"+pimax+"' ecf='"+ecf+"' \
                               esp_nsplinenodes="+esp_nsplinenodes+" esen_mlmin=10;"+ \
                           "cp emllist.fits emllist_"+inst+".fits;"+  
                           "srcdisplay boxlistset=emllist_"+inst+".fits imageset="+name+"_imagesmap.fits sourceradius=0.005 \
                               withregionfile=true regionfile=allSources_"+inst+".reg;"
                           #"fv emllist_"+inst+".fits;" 
                           #"ds9 PNMOS12_image_fullsmap.fits -regions load allSources.reg \
                           #    -cmap bb -scale log -zoom 4 -export 'allSources.jpeg' 300;"
                           , shell=True) 

            t = Table.read(workdir+'/'+'emllist_'+inst+'.fits', format='fits')
            #t = t[list(t.columns)]
            t.write(workdir+'/'+'emllist_'+inst+'.csv', overwrite=True)
            df = pd.read_csv(workdir+'/'+'emllist_'+inst+'.csv').drop_duplicates(['X_IMA'])
            allX, allY = np.array(df.X_IMA), np.array(df.Y_IMA)

            self.otherSourcesIndi[obsID][inst] = [(x, y) for x, y in zip(allX, allY)]

            print('\nAll other Sources for obsID {} found.'.format(obsID))
        #print(self.otherSources)    
        return print('\nFinished finding all the other Sources!')


    ##-- function to get Background Circles --##
    def getBackgroundCircles_indi (self, inst='PN'):
        """
        Calls the `backgroundCircles` function in `findOverlap.py` to automatically 
        find 2 background cicles in the data of the instrument passed.

        Args:
            inst (str): name of the instrument of use.

        :Input: `inst_CCD##_image.fits` files containing the extracted image data of the instrument.

        :Requires: Call to `find_otherSources_indi` function for the dictionary of
                   other Sources.

        Raises:
            Error: if `inst` is not in ['PN', 'MOS1', 'MOS2']

        :Output: x, y coords of the background circles alongwith the radii of these circles.
                Units of radii are pixels or Sky Coord units.
                Updates the ``backgroundLocIndi`` dictionaries with the background location
                parameters.

        NOTES
        -----
        ``small-mode`` MOS checking is invalid here since overlap is not been detected.
        The CCD argument in ``findOverlap`` class is named ``pnCCD``.
        """
        print('\nAutomatically detecting the Background Circles.')
        
        #-- set the environment variables --#
        os.environ['SAS_DIR'] = self.sas_dir
        os.environ['HEADAS'] = self.headas
        os.environ['SAS_CCFPATH'] = self.sas_ccfpath
        
        #-- iterating for all the obsIDs --#      
        for obsID in self.obsIDs:
            print(f'\nDetecting {inst} BkgCircs for obsID {obsID}.')
            workdir = self.workdir+'/'+obsID+'/work'

            self.backgroundLocIndi[obsID] = {}  #--initialize empty dict for saving indi inst bkgCircs.
            
            #-- the source location --#
            instCCD = str(self.sourceCCDs[obsID][inst])
            locParams = self.sourceLoc[obsID]
            srcX, srcY = locParams['x'], locParams['y']
            srcCoords = (srcX, srcY)

            #-- coordinates of other sources --#
            otherSrc = self.otherSourcesIndi[obsID][inst]
            
            #-- find the overlap region and background circles --#
            findOverlapObj = findOverlap(workdir=workdir, obsID=obsID, pnCCD=instCCD, srcCoords=srcCoords, otherSrc=otherSrc, \
                                         srcR=self.srcCircRadius, bkgR=self.bkgCircRadius, srcThreshold=self.dSrcThreshold, \
                                         gap=self.bkgCircGap, saveFig=self.saveFig, showFig=self.showFig)
            bCircles, correctSrc, srcR = findOverlapObj.main(inst=inst)
            
            #-- save the background circle parameters --#
            bLocParams = {}
            for i, bCircle in enumerate(bCircles):
                xKey, yKey, rKey = [key+str(i+1) for key in ['Bx', 'By', 'Br']]
                bLocParams[xKey], bLocParams[yKey], bLocParams[rKey] = bCircle[0], bCircle[1], bCircle[2]
            self.backgroundLocIndi[obsID][inst] = bLocParams

            #-- save the corrected main source coordinates --#
            self.sourceLoc[obsID]['x'], self.sourceLoc[obsID]['y'] = correctSrc[0], correctSrc[1]
            self.sourceLoc[obsID]['r'] = srcR
            
            print('Background for obsID {} found.'.format(obsID))
          
        #print(self.backgroundLoc)
        return print('\nCompleted finding the Background Circles!')


    ##-- write location parameters to the pickle file --##
    def writeCCDcoordsPickle (self):
        """
        Writes the location parameters to `ccd_coords_info.pickle` file.
        If the file is already present in the directory, then the values for the 
        corresponding keys will be updated.

        21st May 2021: Added loading `ccd_coords_info.pickle` file if it is already
                       present in the directory.
        """
        maindir = self.workdir
        print('\nWriting the location parameters to ccd_coords_info.pickle file.\n')

        #-- iterating for all the obsIDs --#      
        for obsID in self.obsIDs:
            print('Saving pickle for obsID {}.'.format(obsID))
            workdir = maindir+'/'+obsID+'/work'
            fname = workdir+'/'+'ccd_coords_info.pickle'

            #-- load the `ccd_coords_info.pickle` if already present --#
            print(self.backgroundLoc[obsID])
            if os.path.isfile(fname):
                result = pickle.load(open(fname, 'rb'))
                    #-- create a backup file --#
                subprocess.run(f"cd '{workdir}';"+ \
                               "cp ccd_coords_info.pickle ccd_coords_info.bak;"
                               , shell=True)
            else:
                result = {}
            
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
            print(self.backgroundLoc[obsID])

            outfile = open(fname, 'wb')
            pickle.dump(result, outfile)
            outfile.close()
            print(result['backgroundLoc'])

        return print('\nWrote location parameters to the pickle file.')


    ##-- update location parameters in the pickle file --##
    def updateCCDcoordsPickle (self, inst=None):
        """
        Writes the background location parameter dictionary for individual instruments
        ``backgroundLoc_indi``  to `ccd_coords_info.pickle` file.
        If the file is already present in the directory, then the values for the 
        corresponding keys are updated.
    
        Args:
            inst (str): name of the instrument of use.

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
                subprocess.run(f"cd '{workdir}';"+ \
                               "cp ccd_coords_info.pickle ccd_coords_info.bak;"
                               , shell=True)
            else:
                self.writeCCDcoordsPickle()
                self.updateCCDcoordsPickle(inst=inst)

            #-- save ``backgroundLoc_indi`` dictionary --#
            if result.get('backgroundLoc_indi', None) is None:
                result['backgroundLoc_indi'] = {}
            result['backgroundLoc_indi'][inst] = self.backgroundLocIndi[obsID][inst]

            #-- save ``otherSources_indi`` dictionary --#
            if result.get('otherSources_indi', None) is None:
                result['otherSources_indi'] = {}
            result['otherSources_indi'][inst] = self.otherSourcesIndi[obsID][inst]

            #-- dump to the pickle file --#
            outfile = open(fname, 'wb')
            pickle.dump(result, outfile)
            outfile.close()

        return print('\nWrote location parameters to the pickle file.')

    
    ##-- read the pickle file for location parameters --##
    def readCCDcoordsPickle (self):
        """
        Reads `ccd_coords_info.pickle` file which contains the location parameters.
        """
        workdir = self.workdir
        print('\nLoading the ccd_coords_info.pickle file.')

        #-- iterating for all the obsIDs --#      
        for obsID in self.obsIDs:
            print('\nReading pickle for obsID {}.'.format(obsID))
            fname = workdir+'/'+obsID+'/work/ccd_coords_info.pickle'

            #-- check for the pickle file --#
            if not os.path.isfile(fname):
                print(f'\nFile {fname} not found!\n')
                self.badObs.append(obsID)
                continue

            #-- load the ccd_coords_info.pickle --#
            result = pickle.load(open(fname, 'rb'))

            #-- get the location and other parameters --#
            self.sourceCCDs[obsID] = result['sourceCCDs']
            self.sourceLoc[obsID] = result['sourceLoc']
            self.backgroundLoc[obsID] = result['backgroundLoc']
            self.otherSources[obsID] = result['otherSources']
            self.smallMode[obsID]  = result['smallMode']
            if result['badObs']:
                self.badObs.append(obsID)
            
        #-- remove badObs from obsIDs list to use --#
        obsIDs, badObs = set(self.obsIDs), set(self.badObs)
        self.obsIDs = list(obsIDs-badObs)

        return print('Read location parameters from the pickle file.')

           
    ##-- function to extract final Source and Background Event lists --##
    def extract_srcBkg_eventLists (self):
        """
        Extracts Source and Background Event list from the combined PNMOS12 data file,
        using the Source and Background location parameters given as Input.
        These Event lists will then be used for getting source images in JPEG using DS9.

        19th May 2021: Added options for Pile-up obsIDs.
        
        :Input: Combined PNMOS12 Event list, Source and Background X, Y in Sky coords and
               rIn, rOut, the radii of the circles.
        :Output: Extracted Source and Background Event lists.
        """
        print('\nInitiating final Source and Background Event lists extraction.')
        
        #-- set the environment variables --#
        os.environ['SAS_DIR'] = self.sas_dir
        os.environ['HEADAS'] = self.headas
        os.environ['SAS_CCFPATH'] = self.sas_ccfpath
        
        #-- iterating for all the obsIDs --#      
        for obsID in self.obsIDs:
            print('\nExtracting Event lists for obsID {}.'.format(obsID))
            workdir = self.workdir+'/'+obsID+'/work'
            
            locParams = self.sourceLoc[obsID]
            srcX, srcY, srcR = str(locParams['x']), str(locParams['y']), str(locParams['r'])
            
            #-- Piled-up cases --#
            srcRin = locParams.get('rIn', None)
            if not self.ignorePileup and srcRin is not None:
                srcRin = str(srcRin)
                srcRout = str(locParams['rOut'])
                print('\nThis observation was piled-up. Using ANNULUS for Source.')
                srcFilterSet = workdir+"/source_PNMOS12_annulus_"+obsID+".evts"
                srcFilterExp = "'((X,Y) in ANNULUS("+srcX+","+srcY+","+srcRin+","+srcRout+"))'"
            else:
                srcFilterSet = workdir+"/source_PNMOS12_"+obsID+".evts"
                srcFilterExp = "'((X,Y) in CIRCLE("+srcX+","+srcY+","+srcR+"))'"
            
            #-- extract Source Event list --#
            subprocess.run(f"cd '{workdir}';"+ \
                           ". $HEADAS/headas-init.sh;"+ \
                           ". $SAS_DIR/setsas.sh;"+ \
                           #"fv PNMOS12.evts;"+ \
                           "evselect table="+workdir+"/PNMOS12.evts"+ \
                               " withfilteredset=yes filteredset="+srcFilterSet+ \
                               " keepfilteroutput=yes expression="+srcFilterExp+";"
                           #"fv source_PNMOS12_"+obsID+".evts;"+ \
                           , shell=True)
            
            bLocParams = self.backgroundLoc[obsID]
            Bx1, By1, Br1 = str(bLocParams['Bx1']), str(bLocParams['By1']), str(bLocParams['Br1'])
            Bx2, By2, Br2 = str(bLocParams['Bx2']), str(bLocParams['By2']), str(bLocParams['Br2'])
                           
            #-- extract the Background Event list --#
            subprocess.run(f"cd '{workdir}';"+ \
                           ". $HEADAS/headas-init.sh;"+ \
                           ". $SAS_DIR/setsas.sh;"+ \
                           "evselect table="+workdir+"/PNMOS12.evts"+ \
                               " withfilteredset=yes filteredset="+workdir+"/background_PNMOS12_"+obsID+".evts"+ \
                               " keepfilteroutput=yes expression='((X,Y) in CIRCLE("+ \
                               Bx1+","+By1+","+Br1+"))||((X,Y) in CIRCLE("+Bx2+","+By2+","+Br2+"))';"
                           #"fv background_PNMOS12_"+obsID+".evts;"
                           , shell=True)
            print('\nExtracting Event lists for obsID {} finished.'.format(obsID))
            
        return print('\nExtraction of final Source and Background Events lists is Successful!')
        

    ##-- function to extract PN source light curve for small mode --##
    def _extract_sMode_srcPNlc (self, obsID):
        """
        Extracts the Source PN light curve for small mode obsIDs.
        The PNMOS12 Source and Background Event lists obtained from 'extract_srcBkg_eventLists'
        have the Source circle in the PNMOS12 overlap region and the Background circles in the 
        PN region for small mode obsIDs. 
        Multiplication of the Background light curve by the Ratio of mean PNMOS12 Source count rate
        and mean PN Source count rate, scales up the Background light curve for small mode obsIDs,
        such that as if they weren't in small mode and were extracted from the overlap region.

        :Input: Source location parameters, PN_CCD##.evts, combined GTI file.
        :Output: Source PN count rate, as an array.
        """
        print('\nExtracting the Source PN Event list and light curve.')
        
        #-- set the environment variables --#
        os.environ['SAS_DIR'] = self.sas_dir
        os.environ['HEADAS'] = self.headas
        os.environ['SAS_CCFPATH'] = self.sas_ccfpath
        
        workdir = self.workdir+'/'+obsID+'/work'
        
        #-- get parameters --#
        pnCCD = str(self.sourceCCDs[obsID]['PN'])
        locParams = self.sourceLoc[obsID]
        srcX, srcY, srcR = str(locParams['x']), str(locParams['y']), str(locParams['r'])
            
        #-- Piled-up cases --#
        srcRin = locParams.get('rIn', None)
        if not self.ignorePileup and srcRin is not None:
            srcRin = str(srcRin)
            srcRout = str(locParams['rOut'])
            srcFilterSet = workdir+"/source_PN_annulus_"+obsID+".evts"
            srcFilterExp = "'((X,Y) in ANNULUS("+srcX+","+srcY+","+srcRin+","+srcRout+"))'"
        else:
            srcFilterSet = workdir+"/source_PN_"+obsID+".evts"
            srcFilterExp = "'((X,Y) in CIRCLE("+srcX+","+srcY+","+srcR+"))'"
            
        #-- extract Source Event list --#
        subprocess.run(f"cd '{workdir}';"+ \
                       ". $HEADAS/headas-init.sh;"+ \
                       ". $SAS_DIR/setsas.sh;"+ \
                       #"fv PN_CCD"+pnCCD+".evts;"+ \
                       "evselect table="+workdir+"/PN_CCD"+pnCCD+".evts"+ \
                           " withfilteredset=yes filteredset="+srcFilterSet+ \
                           " keepfilteroutput=yes expression="+srcFilterExp+";"
                       #"fv source_PN_"+obsID+".evts;"+ \
                       "evselect table=source_PN_"+obsID+".evts withrateset=Y \
                           rateset=source_PN_rate.fits maketimecolumn=Y makeratecolumn=Y \
                           expression='(TIME in gti(combinedGTI_"+obsID+".fits)) && (PI in [300:10000])';"
                       #"dsplot table=source_PN_rate.fits x=TIME y=RATE;"+ \
                       #"fv source_PN_rate.fits;"
                       , shell=True)

        #-- source PN count rate --#
        t = Table.read(workdir+'/'+'source_rate.fits', format='fits', hdu=1)
        t.write(workdir+'/'+'source_rate.csv', overwrite=True)
        df = pd.read_csv(workdir+'/'+'source_rate.csv')

        time, rate, error = df.TIME, df.RATE, df.ERROR

        print('\nFinished small mode Source PN extraction.')
        return (time, rate)


    ##-- function to obtain the light curves --##
    def obtain_lightCurves (self):
        """
        Obtains the light curves from the Source and Background Event lists.
        
        :Input: Source and Background Event lists.
        :Output: Corresponding light curves.

        By default, the background light curve is not binned.
        """
        print('\nFinally obtaining the Source and Background light curves.')

        ##-- the binning function --##
        def __binData (x, y, binSize=25):
            """
            Returns the bins of array y, binned with the binSize between min(x) and max(x).
            Here, array x is Time. 
            The increment in x is in seconds, so the binSize and the difference in interval are same.
            """
            if type(binSize)!=int:
                binSize = int(binSize)

            binned, binTimes = [], []
            intervals = np.arange(min(x), max(x), binSize)

            for i in np.arange(0, len(x), binSize):
                if i==len(x)-binSize:
                    Bins = y[i:] 
                    binned.append( np.mean(Bins) )
                    binTimes.append( np.mean(x[i:]) )
                else:
                    Bins = y[i:i+1]
                    binned.append( np.mean(Bins) )
                    binTimes.append( np.mean(x[i:i+1]) )

            binTimes, binned = np.array(binTimes), np.array(binned)
            return (binTimes, binned)


        #-- set the environment variables --#
        os.environ['SAS_DIR'] = self.sas_dir
        os.environ['HEADAS'] = self.headas
        os.environ['SAS_CCFPATH'] = self.sas_ccfpath

        binSize = self.lcBinSize
        binBkglc = strToBool(self.binBkglc)
        gti_indiThreshold = float(self.gti_indiThreshold)
        
        #-- iterating for all the obsIDs --#      
        for obsID in self.obsIDs:
            print('\nGetting light curves for obsID {}.'.format(obsID))
            workdir = self.workdir+'/'+obsID+'/work'

            smallMode = self.smallMode[obsID]
            locParams, bLocParams = self.sourceLoc[obsID], self.backgroundLoc[obsID]
            #print(locParams, bLocParams)
            
            #-- Piled-up cases --#
            srcRin = locParams.get('rIn', None)
            if not self.ignorePileup and srcRin is not None:
                srcFilterSet = workdir+"/source_PNMOS12_annulus_"+obsID+".evts"
                print('\nThis observation was piled-up. Using ANNULUS for Source.')
            else:
                srcFilterSet = workdir+"/source_PNMOS12_"+obsID+".evts"
            
            #-- extract Source and Background light curve --#
            subprocess.run(f"cd '{workdir}';"+ \
                           ". $HEADAS/headas-init.sh;"+ \
                           ". $SAS_DIR/setsas.sh;"+ \
                           "evselect table="+srcFilterSet+" withrateset=Y \
                               rateset=source_rate.fits maketimecolumn=Y makeratecolumn=Y \
                               expression='(TIME in gti(combinedGTI_"+obsID+".fits)) && (PI in [300:10000])';"
                           #"dsplot table=source_rate.fits x=TIME y=RATE;"+ \
                           #"fv source_rate.fits;"
                           "evselect table=background_PNMOS12_"+obsID+".evts withrateset=Y \
                               rateset=background_rate.fits maketimecolumn=Y makeratecolumn=Y \
                               expression='(TIME in gti(combinedGTI_"+obsID+".fits)) && (PI in [300:10000])';"
                           #"dsplot table=background_rate.fits x=TIME y=RATE;"+ \
                           #"fv background_rate.fits;"
                           , shell=True)

            #-- source rates --#
            tSrc = Table.read(workdir+'/'+'source_rate.fits', format='fits', hdu=1)
            tSrc.write(workdir+'/'+'source_rate.csv', overwrite=True)
            dfSrc = pd.read_csv(workdir+'/'+'source_rate.csv')
            
            #-- background rates --#
            tBkg = Table.read(workdir+'/'+'background_rate.fits', format='fits', hdu=1)
            tBkg.write(workdir+'/'+'background_rate.csv', overwrite=True)
            dfBkg = pd.read_csv(workdir+'/'+'background_rate.csv')

            #-- GTIs --#
            tgti = Table.read(workdir+'/'+'source_rate.fits', format='fits', hdu=3)
            tgti.write(workdir+'/'+'source_GTI.csv', overwrite=True)
            dfgti = pd.read_csv(workdir+'/'+'source_GTI.csv')
            gtiStart, gtiStop = dfgti.START, dfgti.STOP
            
            srcTime, srcRate, srcError = dfSrc.TIME, dfSrc.RATE, dfSrc.ERROR
            bkgTime, bkgRate, bkgError = dfBkg.TIME, dfBkg.RATE, dfBkg.ERROR

            srcMean = np.round(np.mean(srcRate), 5)

            #-- scale up background if small mode --#
            if smallMode:
                print('\nFinding the mean PN only Source count rate, because obsID in small mode.')
                srcPNtime, srcPNrate = self._extract_sMode_srcPNlc(obsID=obsID)
                srcPNmean = np.round(np.mean(srcPNrate), 5)
                factor = srcMean/srcPNmean
                bkgRate = bkgRate * factor
                print('\nMultiplying background rate by (srcPNMOS12/srcPN) count ratio of', np.round(factor, 3))

            #-- find src/bkg area ratio --#
            r, Br1, Br2 = float(locParams['r']), float(bLocParams['Br1']), float(bLocParams['Br2'])
            arRatio = (np.pi*r**2) / (np.pi*Br1**2 + np.pi*Br2**2)
            print('\nSrc-to-Bkg area ratio: ', np.round(arRatio, 3))
            print('Mean Count Rate, Background:{}, Source:{}'.format(np.round(np.mean(bkgRate), 3), srcMean))

            #-- normalize background wrt area ratio --#
            bkgRate = bkgRate*arRatio
            bkgMean = np.round(np.mean(bkgRate), 5)
            print('Mean area normalized Background count rate:\t', bkgMean)

            #-- subtract background from source --#
            srcRate = srcRate - bkgMean
            print('Mean background corrected Source count rate:\t', np.round(np.mean(srcRate), 5))
            
            fig, ax = plt.subplots(1, 1, figsize=(10, 5))
            #ax.plot(srcTime, srcRate, '-', color='green')
            #ax.plot(bkgTime, bkgRate, '-', color='red')

            #-- bin the light curves --#
            if binSize == 'linear':
                srcBinTime, srcBinRate = __binData(srcTime, srcRate, binSize=2*max(srcRate))
                bkgBinTime, bkgBinRate = __binData(bkgTime, bkgRate, binSize=25)
            else:
                binSize = int(binSize)
                srcBinTime, srcBinRate = __binData(srcTime, srcRate, binSize=binSize)
                bkgBinTime, bkgBinRate = __binData(bkgTime, bkgRate, binSize=binSize)
            #ax.plot(srcBinTime, srcBinRate, '-', color='gray')

            #-- iterate for all GTIs --#
            for start, stop in zip(gtiStart, gtiStop):
                    #-- remove short gti --#
                if stop-start < gti_indiThreshold:
                    continue

                    #-- plot the light curve --#
                idxBin = np.logical_and(srcBinTime>=start, srcBinTime<=stop)
                ax.plot(srcBinTime[idxBin], srcBinRate[np.where(idxBin)], '-', color='black')

                if binBkglc:
                    idxBkg = np.logical_and(bkgBinTime>=start, bkgBinTime<=stop)
                    ax.plot(bkgBinTime[idxBkg], bkgBinRate[idxBkg], '-', color='red')
                else:
                    idxBkg = np.logical_and(bkgTime>=start, bkgTime<=stop)
                    ax.plot(bkgTime[idxBkg], bkgRate[idxBkg], '-', color='red')

                    #-- start and stop vertical lines --#
                ax.plot([start, start], [0, max(srcBinRate)], '--', color='lime', linewidth=1.5)
                ax.plot([stop, stop], [0, max(srcBinRate)], '--', color='magenta', linewidth=1.5)
            
            #-- finishing the plot --#
            ax.margins(y=0)
            ax.set_title(str(obsID)+': Source & Background', fontsize=15)
            ax.set_xlabel('Time', fontsize=10)
            ax.set_ylabel('Count rate (counts s$^{-1}$)', fontsize=10)
            ax.legend(['bkg corrected src', 'bkg normalized wrt area'], loc='upper right')
            plotAnal.beautifyPlot([ax], tickNum=6)

            plt.tight_layout()
            if self.saveFig:
                plt.savefig(workdir+'/lightcurve_'+str(obsID)+'.png', dpi=300)
            if self.showFig:
                plt.show()
            plt.close()
                           
            print('\nLight curves for obsID {} obtained.'.format(obsID))
            
        return print('\nFinished getting the Source and Background light curves!')


    ##-- function to save the final results --##
    def save_results (self):
        """
        Does two tasks,
            One, saves a pickle file at each obsID directory containing 
                 sourceCCDs, sourceLoc, backgroundLoc and otherSources for it.
            Two, copies the Source, Background Events lists and other output files 
                 from each obsID directory to a results folder in the main directory.
        """
        print('\nLastly saving results to a pickle file for each obsID.')

        maindir = self.workdir

        #-- make the results directory --#
        if not os.path.isdir(maindir+'/results'):
            subprocess.run(f"cd '{maindir}';"+ \
                           "mkdir results/", shell=True)
        resultdir = maindir+'/results'

        #-- iterating for all the obsIDs --#      
        for obsID in self.obsIDs:
            print('\nSaving results for obsID {}.'.format(obsID))
            workdir = maindir+'/'+obsID+'/work'
            
            #-- copy Event lists and other results --#
            srcFiles = glob.glob(workdir+'/*source*.fits')
            bkgFiles = glob.glob(workdir+'/*background*.fits')
            pngFiles = glob.glob(workdir+'/*overlap*.png') + glob.glob(workdir+'/*circles*.png') \
                        + glob.glob(workdir+'/*lightcurve*.png') + glob.glob(workdir+'/*srcBkg_circs*.png')
            jpegFiles = glob.glob(workdir+'/*.jpeg')
            csvFiles = glob.glob(workdir+'/*.csv')
            pickleFile = glob.glob(workdir+'/*ccd*.pickle')

            files = srcFiles + bkgFiles + pngFiles + jpegFiles + csvFiles + pickleFile
            for file in files:
                fname = os.path.basename(file)         #--get file name from the glob path.
                fname = fname.replace('_'+obsID, '')   #--remove obsID from file name, if it is already there.
                fname = obsID +'_'+ fname              #--add the obsID at the start of file name.

                subprocess.run("cp "+file+f" '{resultdir}'/"+fname, shell=True)

            print(f'Save for obsID {obsID} done.')

        return print('\nSaved the result!')
        

   



#################### End of Program #########################
#############################################################


