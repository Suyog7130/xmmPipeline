
###########################################
###       Using xspec in Python         ###
###########################################

"""
10th April 2021:
---
Getting the spectra fits using PyXspec module.
"""

import subprocess
import matplotlib.pyplot as plt

from xspec import *
from plotAnal import plotAnal


def allSpec (workdir):
    
    pnName = workdir+"/PN_CCD4"
    mos1Name = workdir+"/MOS1_CCD1"
    mos2Name = workdir+"/MOS2_CCD1"
    """
    pnS = Spectrum(pnName+"_spectrum_source.fits")
    pnS.background = pnName+"_spectrum_background.fits"
    pnS.response = pnName+".rmf"
    pnS.response.arf = pnName+".arf"
        
    mos1S = Spectrum(mos1Name+"_spectrum_source.fits")
    mos1S.background = mos1Name+"_spectrum_background.fits"
    mos1S.response = mos1Name+".rmf"
    mos1S.response.arf = mos1Name+".arf"
        
    mos2S = Spectrum(mos2Name+"_spectrum_source.fits")
    mos2S.background = mos2Name+"_spectrum_background.fits"
    mos2S.response = mos2Name+".rmf"
    mos2S.response.arf = mos2Name+".arf"
    """
    pnS = Spectrum(pnName+"_spectrum_grouped.fits")
    mos1S = Spectrum(mos1Name+"_spectrum_grouped.fits")
    mos2S = Spectrum(mos2Name+"_spectrum_grouped.fits")
    
    Plot.xAxis = "KeV"
    
    AllData.ignore("bad")
    pnS.ignore("**-0.3 15.0-**")
    mos1S.ignore("**-0.3 15.0-**")
    mos2S.ignore("**-0.3 15.0-**")
    
    #m1 = Model("tbabs*zashift*(powerlaw+bbody)")
    m1 = Model("zashift*powerlaw")
    
    #print(AllModels.sources)
    #print(m1.componentNames)
    #print(m1.zashift.parameterNames)
    
    #m1.powerlaw.norm = 0.4
    #m1.zashift.Redshift = 2.0
    m1.zashift.Redshift.frozen = False
    
    Xset.abund = "wilm"
    
    Fit.nIterations = 100
    Fit.criticalDelta = 1e-1
    Fit.perform()
    
    Plot.device = "/xs"
    #Plot.xLog = True
    #Plot("model")
    #Plot("data", "model", "residuals")
    Plot("data", "residuals")
    
    return True
    
    

def main (instName, workdir, grouped=True, ax=None):
    
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
        S = Spectrum(fname+"_spectrum_grouped.fits")
    else:
        S = Spectrum(fname+"_spectrum_source.fits")
        S.background = fname+"_spectrum_background.fits"
        S.response = fname+".rmf"
        S.response.arf = fname+".arf"
        #Plot.setRebin(minSig=3, maxBins=1096)  #--doesn't work.
    
    Plot.xAxis = "KeV"
    
    AllData.ignore("bad")
    #AllData.ignore("**-0.3 15.-**")
    S.ignore("**-0.3 15.0-**")   
    m1 = Model("tbabs*zashift*(powerlaw)")
    
    #print(AllModels.sources)
    #print(m1.componentNames)
    #print(m1.zashift.parameterNames)
    
    #m1.powerlaw.norm = 0.4
    #m1.zashift.Redshift = 2.0
    m1.zashift.Redshift.frozen = False
    
    Xset.abund = "wilm"
    
    Fit.nIterations = 100
    Fit.criticalDelta = 1e-1
    Fit.perform()
    
    Plot.device = "/xs"
    #Plot.yLog = True
    #Plot("model")
    Plot("ldata", "residuals")
    
    #-- make the matplotlib plot --#
    Sx, Sy = Plot.x(plotWindow=1), Plot.y(plotWindow=1)
    SxErr, SyErr = Plot.xErr(plotWindow=1), Plot.yErr(plotWindow=1)
    foldedS = Plot.model(plotWindow=1)
    resiX, resiY = Plot.x(plotWindow=2), Plot.y(plotWindow=2)
    resiXerr, resiYerr = Plot.xErr(plotWindow=2), Plot.yErr(plotWindow=2)
    
    if ax is None:
        fig, ax = plt.subplots(2, 1, figsize=(5, 10))
        
    ax[0].errorbar(x=Sx, y=Sy, xerr=SxErr, yerr=SyErr, 
                   marker='.', markersize=3, 
                   ls='none', color='black', linewidth=0.5)
    ax[0].plot(Sx, foldedS, drawstyle='steps-pre')
    ax[1].errorbar(x=resiX, y=resiY, xerr=resiXerr, yerr=resiYerr, 
                   marker='.', markersize=3, 
                   ls='none', color='black', linewidth=0.5)
    ax[1].plot([0, resiX[-1]], [0,0], 'c', color='lightgreen')
    ax[0].set_yscale('log')
    
    if ax is None:
        ax[0].set_title('data and folded model')
        ax[1].set_title('residuals')
    
        #plotAnal.beautifyPlot(ax)
        plt.show()
        plt.close()
    
    return True
    
    



if __name__=="__main__":
    obsID = '0810200701'
    workdir = '/media/suyog/DATA/xmm_obs/'+obsID+'/work'
    instName = 'MOS1_CCD1'
    main(instName=instName, workdir=workdir, grouped=True)
    #allSpec(workdir)
    
    
    
    
    
    
#################### End of Program #########################
#############################################################


