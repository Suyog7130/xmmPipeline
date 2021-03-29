
#####################################################
###           plotAnal functionality              ###
#####################################################

import math
import scipy
import statsmodels.api as sm
import numpy as np
import pandas as pd
import seaborn as sns
import statistics as stt
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import pdfkit as pdf

from decimal import Decimal
from scipy import stats
from skimage.metrics import structural_similarity as ssim


#-- plot analysis class --#
class plotAnal:
 
    def __init__ (self, figure, x, y):
        self.figure = figure
        self.x = x
        self.y = y
        
    def linearFit_plotCal (figure, x, y, xloc=0.05, yloc=0.15, fontsize=15):
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
        line = []
        for i in range(len(x)):
            line.append(float(slope)*x[i] + float(intercept))
            
        if slope<0:
            m = '$-$' + str(round(Decimal(slope), 5))[1:]
        else:
            m = str(round(Decimal(slope), 2))
    
        figure.text(xloc, yloc, 'm = ' + m + ', c = ' +
             str(round(Decimal(intercept),2)),
             transform=figure.transAxes, ha='left', va='center', fontsize=fontsize)
        
    def linearFit (figure, x, y, lineDetails=True, xloc=0.05, yloc=0.10, fontsize=15, color='darkorange'):             
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
        line = []
        for i in range(len(x)):
            line.append(float(slope)*x[i] + float(intercept))
        figure.plot(x, line, label="linear fit", color=color)
    
        if slope<0:
            m = '$-$' + str(round(Decimal(slope), 5))[1:]
        else:
            m = str(round(Decimal(slope), 2))
        
        if lineDetails==True:
             figure.text(xloc, yloc, 'm = ' + m + ', c = ' + str(round(Decimal(intercept),2)), \
                         transform=figure.transAxes, ha='left', va='center', fontsize=fontsize)    
        
    def pearson (figure, x, y, line=False, xloc=0.05, yloc=0.10, decPlace=2, rSquared=False, \
                 fontsize=15, linecolor='darkorange', lineDetails=False):
        pearson = scipy.stats.pearsonr(x, y)
        if pearson[0]<0:
            rPear = '$-$' + str(round(Decimal(pearson[0]), decPlace))[1:]
        else:
            rPear = str(round(Decimal(pearson[0]), decPlace))
        figure.text(xloc, yloc, 'r (pearson) = ' + rPear +
             '  (p='+  str(round(Decimal(pearson[1]),2)) + ')',
             transform=figure.transAxes, ha='left', va='center', fontsize=fontsize)
        if line==True:
            if lineDetails==True:
                plotAnal.linearFit(figure, x, y, xloc=xloc, yloc=yloc-0.12, fontsize=fontsize, color=linecolor)
            else:
                plotAnal.linearFit(figure, x, y, lineDetails=False, fontsize=fontsize, color=linecolor)
        if rSquared==True:
            figure.text(xloc, yloc-0.12, 'r$^2$ value = ' + str( round((pearson[0]**2)*100, 2) ) + '%', \
                        transform=figure.transAxes, ha='left', va='center', fontsize=fontsize)
           
    def spearman (figure, x, y, line=False, xloc=0.05, yloc=0.10, decPlace=2, rSquared=False, \
                 fontsize=15, linecolor='darkorange', lineDetails=False):
        spearman = scipy.stats.spearmanr(x, y)
        if spearman[0]<0:
            rSpear = '$-$' + str(round(Decimal(spearman[0]), decPlace))[1:]
        else:
            rSpear = str(round(Decimal(spearman[0]), decPlace))
        figure.text(xloc, yloc, 'r (spearman) = ' + rSpear +
             '  (p='+  str(round(Decimal(spearman[1]),2)) + ')',
             transform=figure.transAxes, ha='left', va='center', fontsize=fontsize)
        if line==True:
            if lineDetails==True:
                plotAnal.linearFit(figure, x, y, xloc=xloc, yloc=yloc-0.12, fontsize=fontsize, color=linecolor)
            else:
                plotAnal.linearFit(figure, x, y, lineDetails=False, fontsize=fontsize, color=linecolor)
        if rSquared==True:
            figure.text(xloc, yloc-0.12, 'r$^2$ value = ' + str( round((pearson[0]**2)*100, 2) ) + '%', \
                        transform=figure.transAxes, ha='left', va='center', fontsize=fontsize)
        
        
    def pearson_spearman (figure, x, y, line=False, xloc=0.05, yloc=0.10, decPlace=2, rSquared=False, \
                          fontsize=15, linecolor='darkorange', lineDetails=False):
        pearson = scipy.stats.pearsonr(x, y)
        spearman = scipy.stats.spearmanr(x, y)
        rPear = str(round(Decimal(pearson[0]), decPlace))
        rSpear = str(round(Decimal(spearman[0]), decPlace))
        if pearson[0]<0:
            rPear = '$-$' + rPear[1:]
        if spearman[0]<0:
            rSpear = '$-$' + rSpear[1:]
        figure.text(xloc, yloc, 'r (pearson) = ' + rPear + \
            '  (p='+  str(round(Decimal(pearson[1]),2)) + ')', \
            transform=figure.transAxes, ha='left', va='center', fontsize=fontsize)
        figure.text(xloc, yloc-0.06, 'r (spearman) = ' + rSpear + \
            '  (p='+  str(round(Decimal(spearman[1]),2)) + ')', \
            transform=figure.transAxes, ha='left', va='center', fontsize=fontsize)
        if line==True:
            if lineDetails==True:
                plotAnal.linearFit(figure, x, y, xloc=xloc, yloc=yloc-0.12, fontsize=fontsize, color=linecolor)
            else:
                plotAnal.linearFit(figure, x, y, lineDetails=False, fontsize=fontsize, color=linecolor)
        if rSquared==True:
            figure.text(xloc, yloc-0.12, 'r$^2$ value = ' + str( round((pearson[0]**2)*100, 2) ) + '%', \
                        transform=figure.transAxes, ha='left', va='center', fontsize=fontsize)

    def rSquared (figure, y, f):
        """
        This can only be there between the y values and the fitted values f from the fitted line.
        """
        yMean = np.mean(y)
        SStot = np.sum([(yi-yMean)**2 for yi in y])
        SSresd = np.sum([(y[i]-f[i])**2 for i in range(len(y))])
        corr = 1 - SStot/SSresd
        return (corr)

    ##-- PSNR function --##
    def PSNR (figure, x, y, xloc=0.05, yloc=0.10, decPlace=2, fontsize=15):
        """
        Calculates Peak Signal-to-Noise Ratio between x and y and plots its value on the figure.
        """
        MSE = sum( np.square(x-y) )/len(x)
        PSNR = 10*np.log( max(x)**2/MSE )
        PSNR = str(round(PSNR, decPlace))
        figure.text(xloc, yloc, 'PSNR = '+PSNR,
                    transform=figure.transAxes, ha='left', va='center', fontsize=fontsize)
        return
      
        
    ##-- calculate SSIM value --##
    def SSIM_value (x, y, decPlace=2, scikit_value=False):
        """
        Calculates the Structural SIMilarity index between two image data given by x and y.
        """
        if len(x.shape)>1 or len(y.shape)>1:
            x[np.isnan(x)] = 0.0
            y[np.isnan(y)] = 0.0
            x, y = x.flatten(), y.flatten()  #--change to 1D arrays.
            x, y = np.abs(x), np.abs(y)      #--change -ve values.
            
        meanX, meanY = np.mean(x), np.mean(y)
        stdX, stdY = np.std(x), np.std(y)
        varX, varY = np.square(stdX), np.square(stdY)
            
        xDiff, yDiff = x-meanX, y-meanY
        xySum = [x*y for x, y in zip(xDiff, yDiff)]
        covarXY = sum(xySum)/(len(xySum))  #--gotta have len-1 here.
            
        k1, k2 = 0.01, 0.03
        if min(y)!=0.0:
            L = np.log2(max(y)) - np.log2(min(y))
        else:
            L = 255
        c1, c2 = (k1*L)**2, (k2*L)**2
            
        num = (2*meanX*meanY + c1) * (2*covarXY + c2)
        den = (meanX**2 + meanY**2 + c1) * (varX + varY + c2)
        SSIM = num/den
        
        SSIM = round(SSIM, decPlace)
                    
        if scikit_value:
            print(len(x), len(y))
            scikit_ssim = ssim(x, y, data_range=max(y)-min(y))
            return SSIM, scikit_ssim
        else:
            return SSIM
        
    ##-- SSIM function --##
    def SSIM (figure, x, y, xloc=0.05, yloc=0.10, decPlace=2, fontsize=15, scikit_value=False):
        """
        Calculates the Structural SIMilarity index between two image data given by x and y.
        The calculated value is then plotted on the figure passed.
        """
        SSIM = plotAnal.SSIM_value(x, y, decPlace=decPlace)
        
        figure.text(xloc, yloc, 'SSIM = '+str(SSIM),
                    transform=figure.transAxes, ha='left', va='center', fontsize=fontsize)
                    
        if scikit_value:
            scikit_ssim = ssim(x, y, data_range=max(y)-min(y))
            return SSIM, scikit_ssim
        else:
            return SSIM

        
    def beautifyPlot (figures, labelsize=11, lengthMajor=10, tickNum=10, tickDirection='out', \
                      lengthMinor=5, minor=False, grid=False, axisColor=None, yTicks=True, xTicks=True):
        """
        If axes subplots is used, with multiple rows and columns, then the ndarray of axes subobjects 
        can be flattened:
            axes.flatten().tolist()
            
        The axisColor takes in a dictionary axisColor with keys as the spine position and values
        as the colour that has to be inserted. 
        Be wary of the american spelling of 'colour', 'color', used in here.
        NOTE: spines can only be used with figure and not subplots!
        """
        if type(figures)==np.ndarray:
            figures = figures.flatten().tolist()
            
        try:
            for i, figure in enumerate(figures):
                """ plt.MaxNLocator(8) can be used to have set number of ticks """
                #if figure.name!='polar':
                #    figure.xaxis.set_major_locator(tck.AutoLocator())  
                #    figure.yaxis.set_major_locator(plt.MaxNLocator(10))
                
                if xTicks==True:
                    if tickNum!=None:
                        figure.xaxis.set_major_locator(tck.MaxNLocator(tickNum))
                    else:
                        figure.xaxis.set_major_locator(tck.AutoLocator())
                    
                if yTicks==True:
                    if tickNum!=None:
                        figure.yaxis.set_major_locator(tck.MaxNLocator(tickNum))
                    else:
                        figure.yaxis.set_major_locator(tck.AutoLocator())
                        
                figure.tick_params(axis='both', which='major', labelsize=labelsize, length=lengthMajor, direction=tickDirection)
                figure.tick_params(axis='both', which='minor', length=lengthMinor, direction=tickDirection)
                
                if minor==True:
                    figure.xaxis.set_minor_locator(tck.AutoMinorLocator())
                    figure.yaxis.set_minor_locator(tck.AutoMinorLocator())
                    if tickDirection!=None:
                        figure.tick_params(axis='both', which='minor')
                        
                if grid==True:
                    figure.grid(True)

                if axisColor!=None:
                    where = list(axisColor.keys())[i]
                    color = list(axisColor.values())[i]
                    figure.spines[where].set_color(color)
                    figure.tick_params(axis='y', colors=color)
                    if minor==True:
                        figure.tick_params(axis='y', which='minor', colors=color)
                    
        except (TypeError, AttributeError) as e:
            width = len(str(e))+4
            #message = 'beautifyPlot takes a list of axes'.center(width, ' ')
            message = str(e).center(width, ' ')
            print('\n\t\t'+'*'*(width+4))
            print(f'\t\t**{message}**')
            print('\t\t'+'*'*(width+4))



################## End of Program #####################################################
#######################################################################################



