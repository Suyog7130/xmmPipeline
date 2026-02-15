
##############################
###      RGS Pipeline      ###
##############################

"""
Python routine to reduce the RGS data, obtain the Spectra and the light curve from it.
The SAS threads to be used, in order of execution are:
    1. https://www.cosmos.esa.int/web/xmm-newton/sas-thread-rgs
       Main thread to reduce the RGS data. Similar to 'emproc' and 'epproc'.
       Should do all the work and obtain the Spectra and the light curves.
    2. https://www.cosmos.esa.int/web/xmm-newton/sas-thread-rgs2
       To correct the Source coordinates and the extraction masks, if the above 
       step was incorrect.
    3. https://www.cosmos.esa.int/web/xmm-newton/sas-thread-timing
       For extraction of RGS light curves.
"""

import os
import sys
import subprocess
import logging

import glob
import pickle
import argparse
import numpy as np
import pandas as pd

#from xmmPipeline import xmmObj
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


##-- the RGS object class --##
class rgsObj:

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


    def _check_dir(self, path):
        """
        Checks if the given path exists and is writable. If not, it creates the directory and checks permissions.
        """
        if not os.path.isdir(path):
            os.makedirs(path, exist_ok=True)
        if not (os.access(path, os.W_OK) and os.access(path, os.X_OK)):
            logging.error(f"Error: You do not have write/execute permissions in {path}.")
            sys.exit(1)

    ##-- function to find the obsIDs --##
    def findObsIDs (self):
        """
        :Input: Coordinates of the object and path to work directory.
        :Output: List of obsIDs.
        """
        ra, dec, workdir = str(self.ra), str(self.dec), self.workdir
        logging.info(f'\nLooking for obsIDs at RA={ra} and DEC={dec}\nWORKDIR is set at {workdir}')
        
        #-- check if workdir exists --#
        self._check_dir(workdir)
            
        #-- check if browse_extract_wget.pl file exists --#
        if not os.path.isfile(workdir+"/"+"browse_extract_wget.pl"):
            logging.info('\nbrowse_extract_wget.pl not found.')
            subprocess.run("cd "+workdir+";"+
                           "wget -q https://heasarc.gsfc.nasa.gov/FTP/heasarc/software/web_batch/browse_extract_wget.pl", shell=True)
            
            logging.info('browse_extract_wget.pl downloaded.')
            logging.info('\nPlease check the PERL path in the file. If required, correct the path given in first line and save the file.')
            subprocess.run("cd "+workdir+";"+
                           "touch browse_extract_wget.pl", shell=True)
            
        #-- download and save parts of xmmmaster table --##
        subprocess.run("cd "+workdir+";"+ \
                       "chmod +x browse_extract_wget.pl;"+ \
                       "./browse_extract_wget.pl table=xmmmaster position='"+ ra+","+dec+ \
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

        #-- copy and append objName to the file --#
        subprocess.run("cd "+workdir+";"+ \
                       "cp obsIDs_list.dat obsIDs_list_"+objName+".dat;", shell=True) 
        
        return logging.info('\nFound '+str(len(obsIDs))+' obsIDs for the object at given position.')


    ##-- read the pickle file for location parameters --##
    def readPickleFile (self):
        """
        Reads the already obtained pickle file containing location parameters.
        """
        logging.info('\nLoading the Pickle file obtained from xmmPipeline.')

        workdir, objName = self.workdir, self.objName

        #-- load the epicObj --#
        epicObj = pickle.load(open(workdir+'/'+'output_xmmObj_'+objName+'.pickle', 'rb'))

        #-- get the location and other parameters --#
        self.badObs = epicObj.badObs
        self.sourceCCDs = epicObj.sourceCCDs
        self.sourceLoc, self.backgroundLoc = epicObj.sourceLoc, epicObj.backgroundLoc
        self.smallMode, self.otherSources = epicObj.smallMode, epicObj.otherSources
        
        #-- remove badObs from obsIDs list to use --#
        obsIDs, badObs = set(self.obsIDs), set(self.badObs)
        self.obsIDs = list(obsIDs-badObs)

        logging.info('Read location parameters from the xmmPipeline Pickle file.')


    ##-- function to reduce the RGS data --##
    def reduceRGSdata (self):
        """
        Reduces the RGS data using 'rgsproc' to obtain Spectra. 
        A folder names 'RGS/' is created in the work folder of each obsID directory.
        See: https://www.cosmos.esa.int/web/xmm-newton/sas-thread-rgs

        :Input: None.
        :Output: Source and Background RGS Spectra and the associated response matrices.
        """
        logging.info('\nInitializing RGS data reduction.')

        maindir = self.workdir
        
        #-- set the environment variables --#
        os.environ['SAS_DIR'] = self.sas_dir
        os.environ['HEADAS'] = self.headas
        os.environ['SAS_CCFPATH'] = self.sas_ccfpath
        
        #-- iterating for all the obsIDs --#      
        for obsID in self.obsIDs:
            logging.info(f'\nReducing RGS data for obsID {obsID}.')

            #-- set SAS_ODF environment variable --#
            os.environ['SAS_ODF'] = maindir+'/'+obsID+'/ODF'

            #-- check for RGS workdir --#
            workdir = maindir+'/'+obsID+'/work/RGS'
            if not os.path.isdir(workdir):
                subprocess.run("mkdir "+workdir+";", shell=True)

            #-- run rgsproc commands --#
            subprocess.run("cd "+workdir+";"+ \
                           ". $HEADAS/headas-init.sh;"+ \
                           ". $SAS_DIR/setsas.sh;"+ \
                           "cifbuild;"+ \
                           '''export SAS_CCF="`pwd`/ccf.cif";'''+ \
                           #"echo $SAS_CCF;"+ \
                           "odfingest;"+ \
                           '''export SAS_ODF="`pwd`/` ls *SUM.SAS`";'''+ \
                           #"echo $SAS_ODF;"+ \
                           "rgsproc;" 
                           #"fv *R1*SRCLI*;"
                           , shell=True)
            logging.info(f'\nData reduction for obsID {obsID} finished.')
            
        logging.info('\nCompleted RGS data reduction.')



    ##-- function to save the final results --##
    def save_results (self):
        """
        Does two tasks,
            One, saves a pickle file at each obsID directory containing 
                 sourceCCDs, sourceLoc, backgroundLoc and otherSources for it.
            Two, copies the Source, Background Events lists and other output files 
                 from each obsID directory to a results folder in the main directory.
        """
        logging.info('\nLastly saving results to a pickle file for each obsID.')

        maindir = self.workdir

        #-- make the results directory --#
        if not os.path.isdir(maindir+'/results'):
            subprocess.run("cd "+maindir+";"+ \
                           "mkdir results/", shell=True)
        resultdir = maindir+'/results'

        #-- iterating for all the obsIDs --#      
        for obsID in self.obsIDs:
            logging.info(f'\nSaving results for obsID {obsID}.')
            workdir = maindir+'/'+obsID+'/work'
            
            #-- save the CCD and coords info in a pickle file --#
            result = {}
            result['sourceCCDs'] = self.sourceCCDs[obsID]
            result['sourceLoc'] = self.sourceLoc[obsID]
            result['backgroundLoc'] = self.backgroundLoc[obsID]
            result['otherSources'] = self.otherSources[obsID]
            result['smallMode'] = self.smallMode[obsID]

            outfile = open(workdir+'/'+'ccd_coords_info.pickle', 'wb')
            pickle.dump(result, outfile)
            outfile.close()
            
            #-- copy Event lists and other results --#
            srcFiles = glob.glob(workdir+'/*source*.fits')
            bkgFiles = glob.glob(workdir+'/*background*.fits')
            pngFiles = glob.glob(workdir+'/*.png')
            jpegFiles = glob.glob(workdir+'/*.jpeg')
            csvFiles = glob.glob(workdir+'/*.csv')

            files = srcFiles + bkgFiles + pngFiles + jpegFiles + csvFiles
            for file in files:
                fname = os.path.basename(file)         #--get file name from the glob path.
                fname = fname.replace('_'+obsID, '')   #--remove obsID from file name, if it is already there.
                fname = obsID +'_'+ fname              #--add the obsID at the start of file name.

                subprocess.run("cp "+file+" "+resultdir+"/"+fname, shell=True)

            logging.info(f'Save for obsID {obsID} done.')

        logging.info('\nSaved the result!')


##-------------------------------------------------------------------------------------------##


##-- the main function --##
def main (args):
    
    #-- create an object of class spectra --#
    obj = rgsObj(ra=args.ra, dec=args.dec, workdir=args.workdir, \
                 sas_dir=args.sas_dir, headas=args.headas, sas_ccfpath=args.sas_ccfpath)
    
    #-- get obsIDs and the objName --#
    try:
        obj.findObsIDs()
    except KeyError:
        logging.error('\nNo obsIDs found for the given location. Confirm that you are connected to the Internet!')

    #-- check if obsID has been given --#
    if args.obsIDs!=None:
        obj.obsIDs = args.obsIDs
        logging.info('Using the obsIDs passed.')

    #-- run the spectra functions --#
    #obj.readPickleFile()
    obj.reduceRGSdata()

    logging.info('\nHurray! The Method Succesfully ran.')
    #if len(obj.badObs)!=0:
    #    logging.info(f'These obsIDs were excluded from analysis: {obj.badObs}')

    return True


if __name__=="__main__":
    
    description = 'Program to reduce RGS and obtain the Spectra.'
    
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

