
###########################################
###       Using xspec in Python         ###
###########################################

"""
10th April 2021:
---
Getting the spectra fits using PyXspec module.

13th April 2021:
---
Adding ArgParse so that this code can be used by ``spectra.py``
to fit spectra.
Nope. This doesn't actually work.

20th April 2021:
---
It works fine now. See the notes for Google Docs for more info.

08 May 2021:
---
Lots of fitting to the Spectra thing remains.
Note: Gotta make sure that only the available grouped Spectra files are looked for in `pyXspec.py`
      when `all` is passed as the `instName`, since Small-mode obsIDs do not have the MOS Spectra.

20&21 May 2021:
---
Completing the pending work of fitting Models to the Spectra.
Two criterions to be checked for each obsIDs:
    - whether obsID is Piled-up? DONE! 
    - whether obsID is in Small-mode? DONE!
Piled-up cases are already taken care of by `epicSpectra.py` is `ignorePileup` flag is OFF.
For these cases, `spectrum_grouped.fits` will contain `spectrum_source_annulus`.

NOTE: The xspec Model object cannot be pickled or dumped into a JSON file.
So, I think this automation can be done later and for now I can manually do the specModel fitting.

FOUND Out how the output can be accessed, logged, saved and used otherwise!
Bingo!
"""

import argparse
import subprocess
import matplotlib.pyplot as plt

from xspec import *

from epicObj import strToBool, printErrorMessage
from plotAnal import plotAnal


##-- fit Spectra for all instruments together --##
def allSpec (workdir, obsID, model="tbabs*zashift*(powerlaw)", \
             modelParams=None, smallMode=False, grouped=True, \
             saveFig=True, showFig=True):

    #-- starting log --#
    logFile = Xset.openLog('xspeclog.txt')
    logFile = Xset.log
    
    #-- `xspec data` --#
    if grouped:
        pnS = Spectrum(workdir+"/PN_spectrum_grouped_"+obsID+".fits")
        if not smallMode:
            mos1S = Spectrum(workdir+"/MOS1_spectrum_grouped_"+obsID+".fits")
            mos2S = Spectrum(workdir+"/MOS2_spectrum_grouped_"+obsID+".fits")
    else:
        pnS = Spectrum(workdir+"/PN_spectrum_source_"+obsID+".fits")
        pnS.background = workdir+"/PN_spectrum_background_"+obsID+".fits"
        pnS.response = workdir+"/PN_"+obsID+".rmf"
        pnS.response.arf = workdir+"/PN_"+obsID+".arf"
        
        if not smallMode:
            mos1S = Spectrum(workdir+"/MOS1_spectrum_source_"+obsID+".fits")
            mos1S.background = workdir+"/MOS1_spectrum_background_"+obsID+".fits"
            mos1S.response = workdir+"/MOS1_"+obsID+".rmf"
            mos1S.response.arf = workdir+"/MOS1_"+obsID+".arf"
            
            mos2S = Spectrum(workdir+"/MOS2_spectrum_source_"+obsID+".fits")
            mos2S.background = workdir+"/MOS2_spectrum_background_"+obsID+".fits"
            mos2S.response = workdir+"/MOS2_"+obsID+".rmf"
            mos2S.response.arf = workdir+"/MOS2_"+obsID+".arf"
    
    #-- `xspec setplot energy` --#
    Plot.xAxis = "KeV"
    
    #-- `xspec ignore` --#
    AllData.ignore("bad")
    pnS.ignore("**-0.3 1.5-**")
    if not smallMode:
        mos1S.ignore("**-0.3 1.5-**")
        mos2S.ignore("**-0.3 1.5-**")
    
    #-- `xspec model` --#
    #m1 = Model("tbabs*zashift*(powerlaw+bbody)")
    m1 = Model(model)

    #-- input model parameters --#
    """ print(model, modelParams)
    #print(AllModels.sources)
    print(m1.componentNames)
    #print(m1.zashift.parameterNames)
    for cName in m1.componentNames:
        comp = getattr(m1, cName)
        print(comp)
        print(comp.parameterNames) """
    for i in modelParams.keys():
        if type(i) == int:
            printErrorMessage(i)
            param = m1(i)            #--find the `i`th parameter object.
            param.values = modelParams[i]   #--assign value from the modelParams dict.
            #print(param.values)
            freeze = modelParams.get('freeze', None)
            if freeze is not None and i in freeze:
                printErrorMessage(i)
                param.frozen = True

    
    #m1.setPars(modelParams)
    
    """ for modelPart in model.split('*'):
        #if modelPart == 'powerlaw':
        #    m1.powerlaw.norm = 0.4
        if modelPart == "zashift":
            #m1.zashift.Redshift = 2.0
            m1.zashift.Redshift.frozen = False
        if modelPart == "bbody":
            m1.bbody.kT = 0.05 """
    
    #-- `xspec abund` --#
    Xset.abund = "wilm"
    
    #-- `xspec fit` --#
    Fit.nIterations = 100
    Fit.criticalDelta = 1e-1
    Fit.perform()

    #-- `xspec lumin, flux` --#
    AllModels.calcFlux("0.3 1.5 0.02")
    AllModels.calcLumin("0.3 1.5 0.02")
    print(pnS.flux, pnS.lumin)
    if not smallMode:
        print(mos1S.flux, mos1S.lumin)
        print(mos2S.flux, mos2S.lumin)

    print(Fit.statistic, Fit.testStatistic, Fit.dof)
    #print(m1.bbody.kT.values[0], m1.bbody.kT.sigma)
    param1 = m1(1)
    print(param1, param1.values[0])
    
    #-- plotting --#
    Plot.device = "/xs"
    #Plot.xLog = True
    #Plot("model")
    #Plot("data", "model", "residuals")
    #Plot("ldata", "residuals", "background")
    Plot("ldata", "del", "background")
    
    #-- make the matplotlib plot --#
    fig, ax = plt.subplots(3, 1, figsize=(10, 10))

    if smallMode:
        plotGroups, colors, labels = [1], ['black'], ['PN']
    else:
        plotGroups, colors = [1, 2, 3], ['black', 'red', 'green']
        labels = ['PN', 'MOS1', 'MOS2']
    
    for pG, color, label in zip(plotGroups, colors, labels):
        Sx, Sy = Plot.x(plotWindow=1, plotGroup=pG), Plot.y(plotWindow=1, plotGroup=pG)
        SxErr, SyErr = Plot.xErr(plotWindow=1, plotGroup=pG), Plot.yErr(plotWindow=1, plotGroup=pG)
        foldedS = Plot.model(plotWindow=1, plotGroup=pG)
        
        resiX, resiY = Plot.x(plotWindow=2, plotGroup=pG), Plot.y(plotWindow=2, plotGroup=pG)
        resiXerr, resiYerr = Plot.xErr(plotWindow=2, plotGroup=pG), Plot.yErr(plotWindow=2, plotGroup=pG)
    
        Bx, By = Plot.x(plotWindow=3, plotGroup=pG), Plot.y(plotWindow=3, plotGroup=pG)
        BxErr, ByErr = Plot.xErr(plotWindow=3, plotGroup=pG), Plot.yErr(plotWindow=3, plotGroup=pG)
        
        ax[0].errorbar(x=Sx, y=Sy, xerr=SxErr, yerr=SyErr, \
                       marker='.', markersize=3, label=label, \
                       ls='none', color=color, linewidth=0.5)
        ax[0].plot(Sx, foldedS, drawstyle='steps-pre', color=color)
        ax[1].errorbar(x=resiX, y=resiY, xerr=resiXerr, yerr=resiYerr, \
                       marker='.', markersize=3, label=label, \
                       ls='none', color=color, linewidth=0.5)
        ax[1].plot([0, resiX[-1]], [0,0], 'c', color='lightgreen')
        ax[2].errorbar(x=Bx, y=By, xerr=BxErr, yerr=ByErr, \
                       marker='.', markersize=3, label=label, \
                       ls='none', color=color, linewidth=1.0)
        #ax[2].plot(Bx, foldedS, drawstyle='steps-pre', color=color, \
        #           linewidth=0.5, alpha=0.75)
        ax[2].errorbar(x=Sx, y=Sy, xerr=SxErr, yerr=SyErr, \
                       marker='.', markersize=3, label=label+' ldata', \
                       ls='none', color=color, linewidth=0.5, alpha=0.5)
        
    for i in range(3):
        ax[i].set_xscale('log')
        ax[i].legend(loc='upper right')
    ax[0].set_yscale('log')
    ax[2].set_yscale('log')
    
    ax[0].set_title('data and folded model')
    ax[1].set_title('del')
    ax[2].set_title('background')
    
    ax[0].set_ylabel('normalized counts s$^{-1}$ KeV$^{-1}$')
    ax[2].set_ylabel('normalized counts s$^{-1}$ KeV$^{-1}$')
    ax[2].set_xlabel('Energy (KeV)')
    
    plt.suptitle(f'{obsID}\n{model}', x=0.05, y=0.98, horizontalalignment='left')
    plotAnal.beautifyPlot(ax, minor=True, logXformat='scalar', logXminorLabel=True)
    plt.tight_layout()
    if saveFig:
        plt.savefig('EPIC_spectra_'+model.replace('*','-')+'.png', dpi=300)
    if showFig:
        plt.show()
    plt.close()

    #modelFile = open('model.pickle', 'wb')
    #pickle.dump(m1, modelFile)

    Xset.closeLog()   #--close the log.
    
    return True
    
    
##-- fit individual instrument Spectra --##
def indiSpec (instName, workdir, obsID, model="tbabs*zashift*(powerlaw)", grouped=True, \
              saveFig=True, showFig=True):
    
    """
    subprocess.run("cd "+workdir+";"+ \
                   " data l "+instName+"_spectrum_grouped.fits"+ \
                   " response l "+instName+".rmf"+ \
                   " arf l "+instName+".arf"+ \
                   " cpd /xs"+ \
                   " setplot energy"+ \
                   " ignore bad"+ \
                   " ignore **-0.3 15.-**"+ \
                   " plot data;"
                   " model wabs*powerlaw;"
                   " fit 100 1e-1"
                   " setplot rebin 3 4096"
                   " plot data residuals"
                   " error 2.706 1 2 3"
                   , shell=True)
    """
    
    fname = workdir+"/"+instName
    
    if grouped:
        S = Spectrum(fname+"_spectrum_grouped_"+obsID+".fits")
    else:
        S = Spectrum(fname+"_spectrum_source_"+obsID+".fits")
        S.background = fname+"_spectrum_background_"+obsID+".fits"
        S.response = fname+"_"+obsID+".rmf"
        S.response.arf = fname+"_"+obsID+".arf"
        #Plot.setRebin(minSig=3, maxBins=1096)  #--doesn't work.
    
    Plot.xAxis = "KeV"
    
    AllData.ignore("bad")
    #AllData.ignore("**-0.3 15.-**")
    S.ignore("**-0.3 10.0-**")   
    m1 = Model(model)
    
    #print(AllModels.sources)
    #print(m1.componentNames)
    #print(m1.zashift.parameterNames)
    
    for modelPart in model.split('*'):
        #if modelPart == 'powerlaw':
        #    m1.powerlaw.norm = 0.4
        if modelPart == "zashift":
            #m1.zashift.Redshift = 2.0
            m1.zashift.Redshift.frozen = False
        if modelPart == "bbody":
            m1.bbody.kT = 0.05
    
    Xset.abund = "wilm"
    
    Fit.nIterations = 100
    Fit.criticalDelta = 1e-1
    Fit.perform()
    
    Plot.device = "/xs"
    #Plot.yLog = True
    #Plot("model")
    Plot("ldata", "residuals", "background")
    
    #-- make the matplotlib plot --#
    Sx, Sy = Plot.x(plotWindow=1), Plot.y(plotWindow=1)
    SxErr, SyErr = Plot.xErr(plotWindow=1), Plot.yErr(plotWindow=1)
    foldedS = Plot.model(plotWindow=1)
    
    resiX, resiY = Plot.x(plotWindow=2), Plot.y(plotWindow=2)
    resiXerr, resiYerr = Plot.xErr(plotWindow=2), Plot.yErr(plotWindow=2)
    
    Bx, By = Plot.x(plotWindow=3), Plot.y(plotWindow=3)
    BxErr, ByErr = Plot.xErr(plotWindow=3), Plot.yErr(plotWindow=3)
    
    #if ax is None:
        #print('printing the matplotlib plot.')
    fig, ax = plt.subplots(2, 1, figsize=(10, 10))
        
    ax[0].errorbar(x=Sx, y=Sy, xerr=SxErr, yerr=SyErr, \
                   marker='.', markersize=3, \
                   ls='none', color='black', linewidth=0.5)
    ax[0].plot(Sx, foldedS, drawstyle='steps-pre', color='black')
    ax[1].errorbar(x=resiX, y=resiY, xerr=resiXerr, yerr=resiYerr, \
                   marker='.', markersize=3, \
                   ls='none', color='black', linewidth=0.5)
    ax[1].plot([0, resiX[-1]], [0,0], 'c', color='lightgreen')
    ax[0].errorbar(x=Bx, y=By, xerr=BxErr, yerr=ByErr, \
                   marker='.', markersize=3, \
                   ls='none', color='blue', alpha=0.75, linewidth=0.5)
    ax[0].set_xscale('log')
    ax[1].set_xscale('log')
    ax[0].set_yscale('log')
    
    #if ax is None:
    ax[0].set_title('data, background and folded model')
    ax[1].set_title('residuals')
    
    ax[0].set_ylabel('normalized counts s$^{-1}$ KeV$^{-1}$')
    ax[1].set_xlabel('Energy (KeV)')
    
    plt.suptitle(f'{obsID}\n{instName}\n{model}', x=0.05, y=0.98, horizontalalignment='left')
    plotAnal.beautifyPlot(ax, minor=True, logXformat='scalar', logXminorLabel=True)
    plt.tight_layout()
    if saveFig:
        plt.savefig(instName.split('_')[0]+'_spectra_'+model.replace('*','-')+'.png', dpi=300)
    if showFig:
        plt.show()
    plt.close()
    
    return True
    
    



if __name__=="__main__":

    description = 'Program to fit a model to the Spectra.'
    
    parser = argparse.ArgumentParser(description=description)   
    
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
    
    #-- general arguments --#
    parser.add_argument('--workdir', action='store', type=str, default='/media/suyog/DATA/xmm_obs', \
                        help='directory where obsid folders will be stored. (default:%(default)s)')
    parser.add_argument('--obsID', action='store', default=None, #['0810200701'], \
                        help='obsID to fit the Spectra for. (default:%(default)s)')
    parser.add_argument('--instName', action='store', default='PN', \
                        help='the instrument to use. (default:%(default)s)')
    parser.add_argument('--model', action='store', default="tbabs*zashift*(powerlaw)", \
                        help='what model to use for fitting. (default:%(default)s)')
    parser.add_argument('--modelParams', action='store', type=eval, \
                        help='give a dictionary of model parameter values to use.')

    parser.add_argument('--grouped', action='store', default='yes', \
                        help='whether to use grouped spectra data or not? (default:%(default)s)')
    parser.add_argument('--smallMode', action='store', default='no', \
                        help='is the obsID in Small-mode? (default:%(default)s)')

    parser.add_argument('--noSaveFig', action='store_false', default=True, \
                        help='do not save the matplotlib plots? (default:%(default)s)')
    parser.add_argument('--showFig', action='store_true', default=False, \
                        help='show the matplotlib plots or not? (default:%(default)s)')

    #-- parse the arguments --#
    args = parser.parse_args()
    
    obsID = args.obsID
    workdir = '/media/suyog/DATA/xmm_obs/'+obsID+'/work'
    instName, model = args.instName, args.model
    grouped, smallMode = strToBool(args.grouped), strToBool(args.smallMode)
    modelParams = args.modelParams #eval(args.modelParams)
    
    if instName == 'all':
        allSpec(workdir, obsID=obsID, model=model, modelParams=modelParams,
                grouped=grouped, smallMode=smallMode, \
                saveFig=args.noSaveFig, showFig=args.showFig)
    else:
        indiSpec(instName=instName, workdir=workdir, obsID=obsID, model=model, \
                 grouped=grouped, saveFig=args.noSaveFig, showFig=args.showFig)
    
    
#################### End of Program #########################
#############################################################




