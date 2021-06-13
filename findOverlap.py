
###############################
###      Find Overlap       ###
###############################

"""
08-12 February 2021:
---
New Python function to obtain the overlap region.

22 February 2021:
---
Completing the modifications to the codes that were discussed on 13th of this month.

26 February 2021:
---
Had the call with Dheeraj@MIT yester. Finalizing the things about this pipeline now.

10 March 2021:
---
The radius of the Source and Background circles given as input to evselect have to be in Sky Coords 
or Pixel units and not arcsec, as I had been giving. Thus, the returned values from findOverlap.py 
have to be altered.

30 March 2021:
---
Was thinking of using numba to accelerate the code, but since most of the functions 
use matplotlib it is perhaps unfeasible.

29th April 2021:
---
The Background circles were too small in some cases.
Modifying the criterion for background circles radius and threshold.

Instead of taking some random point from with the points that are a threshold distance
away from other Sources, I can take the one with max radii. Currently, the Bkg radius is defined a-priori.

Tried using this, but there are so many complications.
So, have simply added another function to find the max bkg radii possible within the threshold,
to be used after the bkg coords have been found.

I think, I can have Bx1,By1 be found together with the xUseN,yUseN and continue having `__maxBkgRadius`
function for the Radii.
Yeah! This works. However, whether or not optimum radii are obtained, depends on which of the circles
is taken as first and which as the second.

I can also try to have the `xmmObj` class as a meta-class for the `findOverlap` class.
Nope! It's better not to do this because I run `findOverlap.py` for each obsID separately.

30th April 2021:
---
Made the second background point non-random as well.
The program takes a little longer now, but both the circles are large sized.
Gonna check if there's any errors that come up because of this.
Also lowered the Bkg Circle Gap to 2.5 pixels instead of 5 pixles.

13th June 2021:
---
Writing the functions for obtaining individual background circles.
"""

import os
import time
import glob
import argparse
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from astropy.io import fits

from skimage import feature  #--for canny edge filter.
from skimage import measure  #--for find_contours.
from skimage.transform import probabilistic_hough_line
from skimage.draw import polygon2mask, polygon

import cv2 as opencv



def imagePN (workdir, pnCCD):
    subprocess.run("cd "+workdir+";"+ \
                   ". $HEADAS/headas-init.sh;"+ \
                   ". $SAS_DIR/setsas.sh;"+ \
                   "evselect table=PN_CCD"+pnCCD+".evts:EVENTS imagebinning='binSize' \
                       imageset='PN_CCD"+pnCCD+"_image.fits' withimageset=yes \
                       xcolumn='X' ycolumn='Y' ximagebinsize=80 yimagebinsize=80 \
                       expression='#XMMEA_EP && (PI in [300:10000])&&(PATTERN in [0:12])';"
                   , shell=True)
                   
def imageMOS12 (workdir):
    subprocess.run("cd "+workdir+";"+ \
                   ". $HEADAS/headas-init.sh;"+ \
                   ". $SAS_DIR/setsas.sh;"+ \
                   "evselect table=MOS12.evts:EVENTS imagebinning='binSize' \
                       imageset='MOS12_image.fits' withimageset=yes \
                       xcolumn='X' ycolumn='Y' ximagebinsize=80 yimagebinsize=80 \
                       expression='#XMMEA_EM && (PI in [300:10000]) && (PATTERN in [0:12])';"
                   , shell=True)
                   
def imagePNMOS12 (workdir):
    subprocess.run("cd "+workdir+";"+ \
                   ". $HEADAS/headas-init.sh;"+ \
                   ". $SAS_DIR/setsas.sh;"+ \
                   "evselect table=PNMOS12.evts:EVENTS imagebinning='binSize' \
                       imageset='PNMOS12_image_full.fits' withimageset=yes \
                       xcolumn='X' ycolumn='Y' ximagebinsize=80 yimagebinsize=80 \
                       expression='(PI in [300:10000]) && (PATTERN in [0:12])';"
                   , shell=True)
                   
def imagePN_soft (workdir, pnCCD):
    subprocess.run("cd "+workdir+";"+ \
                   ". $HEADAS/headas-init.sh;"+ \
                   ". $SAS_DIR/setsas.sh;"+ \
                   "evselect table=PN_CCD"+pnCCD+".evts:EVENTS imagebinning='binSize' \
                       imageset='PN_CCD"+pnCCD+"_image_soft.fits' withimageset=yes \
                       xcolumn='X' ycolumn='Y' ximagebinsize=20 yimagebinsize=20 \
                       expression='#XMMEA_EP && (PI in [300:1000]) && (PATTERN==0)';"
                   , shell=True)
                   
def imageMOS12_soft (workdir):
    subprocess.run("cd "+workdir+";"+ \
                   ". $HEADAS/headas-init.sh;"+ \
                   ". $SAS_DIR/setsas.sh;"+ \
                   "evselect table=MOS12.evts:EVENTS imagebinning='binSize' \
                       imageset='MOS12_image_soft.fits' withimageset=yes \
                       xcolumn='X' ycolumn='Y' ximagebinsize=20 yimagebinsize=20 \
                       expression='#XMMEA_EM && (PI in [300:1000]) && (PATTERN==0)';"
                   , shell=True)

def imagePNMOS12_soft (workdir):
    subprocess.run("cd "+workdir+";"+ \
                   ". $HEADAS/headas-init.sh;"+ \
                   ". $SAS_DIR/setsas.sh;"+ \
                   "evselect table=PNMOS12.evts:EVENTS imagebinning='binSize' \
                       imageset='PNMOS12_image_soft.fits' withimageset=yes \
                       xcolumn='X' ycolumn='Y' ximagebinsize=20 yimagebinsize=20 \
                       expression='(PI in [300:1000]) && (PATTERN==0)';"
                   , shell=True)


##-- findOverlap Class --##
class findOverlap:

    def __init__ (self, workdir, obsID='NA', pnCCD='4', srcCoords=None, otherSrc=None, gap=2.5, \
                  srcR=None, bkgR=None, srcThreshold=None, saveFig=True, showFig=True):

        self.workdir = workdir
        self.obsID = obsID
        self.pnCCD = pnCCD

        self.overlap = None
        self.imgPNMOS = None
        self.imgMOS = None
        self.imgPN = None
        self.header = None

        self.srcCoords = srcCoords        #--main source coordinates.
        self.otherSrc = otherSrc          #--coordinates of other sources.
        self.srcR = srcR                  #--source circle radius passed by the user or that calculated, arcsec.
        self.bkgR = bkgR                  #--background circle radii passed by the user, arcsec.
        self.srcThreshold = srcThreshold  #--threshold distance from Source passed by the user, arcsec.
        self.gap = gap                    #--gap to have around the circles, in pixels.

        self.bCircle1 = None     #--first found background circle, x, y and radius.
        self.bCircle2 = None     #--second found background circle, x, y and radius.
        self.correctSrc = None   #--corrected main source coordinates.

        self.axlims = {'xlim':None, 'ylim':None}    #--four corners of the PNMOS12 image, for setting the axis limits.

        self.isSmallMode = False
        self.saveFig = saveFig
        self.showFig = showFig


    ##-- function to open the image --##
    def __openImg (self, fname, header=False, printHeader=False):
        """
        Open the required FITS image file.
        Return header as well, if keyword 'header' is set.
        """
        hdu = fits.open(fname)
        img = hdu[0].data
        img = np.log(1000*img+1)
        img = img/np.log(img)
        
        if printHeader:
            for key in header.keys():
                print('{}\t{}\t{}'.format(key, header.comments[key], header[key]))
        if header:
            header = hdu[0].header
            return (img, header)
        return img
        

    ##-- to detect individual collection regions --##    
    def _detect_individual (self, img, wSize=8, percent=0.20, usingCount=False, ax=None, mesh=True):
    
        xDim, yDim = img.shape[0], img.shape[1]
        centerX, centerY = int(xDim/2), int(yDim/2)  
        xlim1, xlim2 = centerX-int(centerX/2), centerX+int(centerX/2)
        ylim1, ylim2 = centerY-int(centerY/2), centerY+int(centerY/2)
    
        wSize = wSize 
        if wSize%2!=0:
            wSize += 1
        frac = int(wSize/2)
    
        percent = percent
        if usingCount==True:
            threshold = percent * (wSize**2)
        else:
            threshold = percent * (wSize**2) * max(img.flatten())
    
        minXpt, maxXpt, minYpt, maxYpt = [xlim2,0], [xlim1,0], [0,ylim2], [0,ylim1]
        for x in range(xlim1, xlim2, wSize):
            for y in range(ylim1, ylim2, wSize):
                if (x<xlim1+frac or x>xlim2-frac) or (y<ylim1+frac or y>ylim2-frac):
                    continue
                window = img[x-frac:x+frac, y-frac:y+frac]
                xl = [x-frac, x-frac, x+frac, x+frac, x-frac]
                yl = [y-frac, y+frac, y+frac, y-frac, y-frac]
                if ax!=None and mesh:
                    ax.plot(xl, yl, '-', color='grey')
            
                if usingCount==True:
                    window = window.flatten()
                    intensity = len(window[window>0.0])
                else:
                    intensity = sum(window.flatten())
            
                if intensity>=threshold:
                    if x<minXpt[0]:
                        minXpt = [x, y]
                    if x>maxXpt[0]:
                        maxXpt = [x, y]
                    if y<minYpt[1]:
                        minYpt = [x, y]
                    if y>maxYpt[1]:
                        maxYpt = [x, y]
                    if ax!=None:
                        ax.fill(yl, xl, fill=True, hatch=3*'//', alpha=0.5, color='brown')

        xA1, yA1 = minXpt[0], minXpt[1]
        xA2, yA2 = maxXpt[0], maxXpt[1]
        xB1, yB1 = minYpt[0], minYpt[1]
        xB2, yB2 = maxYpt[0], maxYpt[1]
    
        xi = [xA1, xB2, xA2, xB1]
        yi = [yA1, yB2, yA2, yB1]
    
        if ax!=None:
            ax.set_xlim([xlim1, xlim2])
            ax.set_ylim([ylim1, ylim2])
            for xCorner, yCorner in zip(xi, yi):
                #print(yCorner, xCorner)
                ax.plot(yCorner, xCorner, color='gray', markersize=10, marker='o')
        #print('\n')
        return (xi, yi)   


    ##-- function to obtain the overlap --##
    def obtainOverlap (self, axes):
        """
        :Input: The workdir where the PN, MOS and PNMOS12 FITS files are located.
        :Output: The detected overlap region if it is larger in size than then MOS region,
                or else the PNMOS12 image data.
        """
        workdir, pnCCD = self.workdir, self.pnCCD
        if not os.path.isfile(workdir+'/'+'PN_CCD'+pnCCD+'_image.fits'):
            imagePN(workdir, pnCCD)
        if not os.path.isfile(workdir+'/'+'MOS12_image.fits'):
            imageMOS12(workdir)
        if not os.path.isfile(workdir+'/'+'PNMOS12_image_full.fits'):
            imagePNMOS12(workdir)
        
        #-- open the images --#
        imgMOS = self.__openImg(workdir+'/'+'MOS12_image.fits')
        imgPN = self.__openImg(workdir+'/'+'PN_CCD'+pnCCD+'_image.fits')
        imgPNMOS, header = self.__openImg(workdir+'/'+'PNMOS12_image_full.fits', header=True)
    
        #-- detect individual collection regions --#
        xiMOS, yiMOS = self._detect_individual(imgMOS, wSize=4, ax=axes[0])
        xiPN, yiPN = self._detect_individual(imgPN, wSize=4, ax=axes[1])
        xiPNMOS, yiPNMOS = self._detect_individual(imgPNMOS, ax=None)
    
        #-- use masked array to find overlap --#
        mask = np.zeros(imgPNMOS.shape)
        xxMOS, yyMOS = polygon(yiMOS, xiMOS)
        xxPN, yyPN = polygon(yiPN, xiPN)
    
        mask[xxMOS, yyMOS] = 0.5
        mask[xxPN, yyPN] = mask[xxPN, yyPN] + 0.5
        mask[mask==0.5] = 0.0
    
        overlap = imgPNMOS*mask.T

        #-- plot the overlap results --#
        plot = axes[0].imshow(imgMOS, cmap='afmhot', origin='lower')
        axes[1].imshow(imgPN, cmap='afmhot', origin='lower')
    
        axes[2].imshow(imgPNMOS, cmap='afmhot', origin='lower')
        for xCorner, yCorner in zip(xiMOS, yiMOS):
            axes[2].plot(yCorner, xCorner, color='gray', markersize=10, marker='o')
        for xCorner, yCorner in zip(xiPN, yiPN):
            axes[2].plot(yCorner, xCorner, color='gray', markersize=10, marker='o')
        #for xCorner, yCorner in zip(xiPNMOS, yiPNMOS):
        #    axes[2].plot(yCorner, xCorner, color='gray', markersize=20, marker='o')

        axes[3].imshow(overlap, cmap='afmhot', origin='lower')
    
            #-- zoom the PNMOS12 plot --#
        axes[2].set_xlim([min(yiPNMOS)-10, max(yiPNMOS)+10])
        axes[2].set_ylim([min(xiPNMOS)-10, max(xiPNMOS)+10])
    
        axes[0].set_title('MOS12', fontsize=12)
        axes[1].set_title('PN', fontsize=12)
        axes[2].set_title('PNMOS12', fontsize=12)
        axes[3].set_title('overlap', fontsize=12)

        #-- pass overlap to other functions --#
        self.overlap, self.imgPNMOS, self.imgMOS, self.imgPN = overlap, imgPNMOS, imgMOS, imgPN
        self.header = header
    
        return True
    
    
    ##-- function to find some background circles --##
    def backgroundCircles (self, axes):
        """
        :Input: The output from 'obtainOverlap()' function aka detected overlap region.,
               and the Source coordinates.
        :Output: Coordinates and Radii of the Background Circles.
        """
        workdir, obsID = self.workdir, self.obsID
        overlap, imgPNMOS, imgMOS, imgPN = self.overlap, self.imgPNMOS, self.imgMOS, self.imgPN
        header = self.header
        srcCoords, otherSrc = self.srcCoords, np.array(self.otherSrc)

        #-- check if overlap region is small --#
        overlap1, imgMOSa, imgPNa = np.nan_to_num(overlap), np.nan_to_num(imgMOS), np.nan_to_num(imgPN)   #--remove invalid values.
        if np.abs( len(overlap1[overlap1>0.0].flatten()) - len(imgMOS[imgMOS>0.0].flatten()) ) < 500: 
            overlap = imgPNMOS
            self.isSmallMode = True
            print('\nNote: obsID in Small Mode.')
    
        ##-- distance between two points --##
        def __euclideanDist (P1, P2):
            """
            Calculates the Euclidean distance between two cartesian points.
        
            :Input: Coordinates of the two points in form of tuples.
            :Output: The distance between the points, float.
        
            Note that the Source coordinates are reversed, so that has to be taken 
            into account when passing the values.
            """
            x1, y1 = P1
            x2, y2 = P2
            return np.sqrt( (x2-x1)**2 + (y2-y1)**2 )

        ##-- distance of a point from a line --##
        def __pDistToLine (P0, line, ax=None):
            """
            Calculates the perpendicular distance of point from a line.
        
            :Input: Two tuples P0 and line. P0 contains the x and y coordinates of the point and
                   line contains the coordinates of the start and end points of the line.
            :Output: The perpendicular distance of the point from the line, float.
        
            First the intersection point Pi of the perpendicular from P0 to the line is calculated.
            The output is then the EuclideanDist between P0 and Pi.
        
            Note that the corner points of the overlap are also reversed in coordinates.
            So, the start and end points of the line have reversed x and y.
            """
            x0, y0 = P0
            x1, y1, x2, y2 = line
        
            m = (y2-y1)/(x2-x1)
            mP = -1/m
        
            xI = ((y1-y0) + mP*x0 - m*x1) / (mP-m)
            yI = y0 + mP * (xI-x0)
        
            if ax!=None:
                ax.plot([x1, x2], [y1, y2], '-', color='blue')
                ax.plot(xI, yI, 'p', markersize=10, color='red')
                ax.plot([x0, xI], [y0, yI], '-', color='violet')
            
            return __euclideanDist(P0, (xI, yI))
        
    
        #-- find corners of the overlap --#
        xiOverlap, yiOverlap = self._detect_individual(overlap, wSize=4, percent=0.05, ax=axes[3], mesh=False)
    
        #-- coordinates of all the points in the overlap --#
        xxOverlap, yyOverlap = polygon(yiOverlap, xiOverlap)
        #axes[3].plot(xxOverlap, yyOverlap, '.', color='pink')
       
        #-- lines made by the overlap corners --#
        cornerLines = []
        for i in range(4):
            if i==3:
                cornerLines.append( (yiOverlap[i], xiOverlap[i], yiOverlap[0], xiOverlap[0]) )
            else:
                cornerLines.append( (yiOverlap[i], xiOverlap[i], yiOverlap[i+1], xiOverlap[i+1]) )
    
        #-- the source coordinates --#
        if srcCoords!=None:
            srcCoords = np.array(srcCoords)/header['CDELT2L']
            xC, yC = srcCoords
            """
            if self.isSmallMode:
                xC, yC = srcCoords
            else:
                yC, xC = srcCoords   #--the axes are opposite here.
            """
            print('\nThe Source coordinate value in log scale pixel is: {}, {}'.format(xC, yC))
        else:
            xC, yC = int(imgPNMOS.shape[0]/2), int(imgPNMOS.shape[1]/2)
        #axes[3].plot(xC, yC, '*', markersize=20, color='white', markeredgecolor='black', markeredgewidth=0.2)

        #-- correct the position of main Source --#
        srcRepeat = np.repeat(np.array([(xC, yC)]), len(otherSrc), axis=0)
        oSrcDist = np.array( list(map(__euclideanDist, otherSrc, srcRepeat)) )
        correctSrc = otherSrc[np.where(oSrcDist==min(oSrcDist))]
        xC, yC = correctSrc[0][0], correctSrc[0][1]
        self.correctSrc = np.array([xC, yC])*header['CDELT2L']
        axes[3].plot(xC, yC, '*', markersize=20, color='white', markeredgecolor='black', markeredgewidth=0.2)

        #-- plot other sources --#
        axes[3].plot(otherSrc[:, 0], otherSrc[:, 1], 'p', markersize=10, color='white', markeredgecolor='black', markeredgewidth=0.2)
    
        #-- threshold distance from the Source --#
        dSrcThreshold = 70  #--value in arcsec.
        if self.srcThreshold != None:
            srcThreshold = float(self.srcThreshold)
            if srcThreshold <= dSrcThreshold:
                dSrcThreshold = srcThreshold
            else:
                print(f'\nProvided threshold value exceed {dSrcThreshold} arcsec. Defaulting to this value.')

        DSource = dSrcThreshold*(1/3600)/header['CDELT2']
        print('\nThe distance between random Background circle center and the Source is set at {} arcsec.'.format(dSrcThreshold))

        #-- radius of the background circles --#
        Br = np.array([dSrcThreshold/2, dSrcThreshold/2 - 10])  #--value in arcsec.
        Br1, Br2 = Br*(1/3600)/header['CDELT2']
        if self.bkgR!=None:
            bkgR = np.array(self.bkgR).astype(float)
            if sum(bkgR) <= sum(Br):
                Br = bkgR
                Br1, Br2 = bkgR*(1/3600)/header['CDELT2']
            else:
                print('\nProvided background circle radii are more than half of threshold distance.')
        #print(Br1, Br2)
    
        #-- dSrcThreshold gives the Source circle radius --#
        srcR = dSrcThreshold/2
        srcToLines = [__pDistToLine((xC, yC), cornerLine) for cornerLine in cornerLines]  #--output is in pixels.
        if min(srcToLines) < srcR*(1/3600)/header['CDELT2']:
            print('\nNote: The Source is very near the overlap edge. Source circle area maybe small.')
            srcR = min(srcToLines)*3600*header['CDELT2']  #--convert to arcsec

        if self.srcR != None:
            srcR_given = float(self.srcR)   #--convert from arcsec.
            if srcR_given <= srcR:
                srcR = srcR_given
            else:
                print(f'\nsrcCircRadius value of {np.round(srcR_given, 2)} arcsec provided exceed the distance to the nearest overlap border. Defaulting to this distance.')

        print(f'The Source radius is set at {np.round(srcR, 2)} arcsec.')
        self.srcR = (srcR*(1/3600)/header['CDELT2'] )*header['CDELT2L']  #--convert pixels to log pixels.

        #-- shortlisting random points --#
        xUse, yUse, xUseN, yUseN = [], [], [], []
        minDistToOtherSrc = Br1
        Bx1, By1 = None, None
        for ptX, ptY in zip(xxOverlap, yyOverlap):
            ptRepeat0 = np.repeat(np.array([(ptX, ptY)]), len(cornerLines), axis=0)
            pDistList = list(map(__pDistToLine, ptRepeat0, cornerLines))
            pDistList.sort()

            #-- check distance from main Source and cornerLines --#
            if __euclideanDist((ptX, ptY), (xC, yC)) > (DSource + self.gap) and pDistList[0] > (Br1 + self.gap):
                    xUse.append(ptX)
                    yUse.append(ptY)

                    #-- check distance from other sources --#
                    ptRepeat = np.repeat(np.array([(ptX, ptY)]), len(otherSrc), axis=0)
                    distList = list(map(__euclideanDist, ptRepeat, otherSrc))
                    distList.sort()
                    if distList[0] > DSource:
                        xUseN.append(ptX)
                        yUseN.append(ptY)
                        if distList[0] >= (minDistToOtherSrc + self.gap):
                            Bx1, By1 = ptX, ptY
 
        #print(len(xUseN), len(yUseN))
        
        ##-- find random background points --##
        def __backgroundPt (xToUse=xUseN, yToUse=yUseN, bkgPt=None, ax=None):
            """ 
            Internal function to select random constrained points.
            At first, the arrays (xUseN, yUseN) are used to find the random points.

            If the number of points in these is small and two random points atleast the 
            threshold distance away from eachother are not found, then the arrays (xUse, yUse)
            are used. 
            """
           
            message = f'\nNOTE:There are too many other sources for obsID {obsID}, such that at the threshold given, no background circles are possible.' \
                        +'\nKindly redo the product extraction for this one later with some changed parameters.' \
                        +'\nFor now the other sources are not taken into account for obtaining the background circles.\n'
            
            #-- check if the arrays are empty --#
            if len(xToUse)<2:
                print(message)
                return __backgroundPt(xToUse=xUse, yToUse=yUse, ax=ax)

            #-- first random point --#
            if bkgPt is None:
                idx1 = np.random.choice( range(len(xToUse)) )
                Bx1, By1 = xToUse[idx1], yToUse[idx1]
            else:
                Bx1, By1 = bkgPt

            #-- second background point --#
            totalBr = Br1 + Br2
            Bx2, By2 = None, None
            for ptX, ptY in zip(xToUse, yToUse):
                dist = __euclideanDist((Bx1, By1), (ptX, ptY))
                if dist >= totalBr:
                    totalBr = dist
                    Bx2, By2 = ptX, ptY

            if Bx2 is not None:
                axes[3].plot(xToUse, yToUse, '.', color='cyan', alpha=0.25)
                return [Bx1, By1, Bx2, By2]
                
            #-- when other sources are too many --#
            else:
                print(message)
                return __backgroundPt(xToUse=xUse, yToUse=yUse, ax=ax)

        if Bx1 is not None:
            Bx1, By1, Bx2, By2 = __backgroundPt(bkgPt=(Bx1, By1), ax=axes[3])
        else:
            Bx1, By1, Bx2, By2 = __backgroundPt(ax=axes[3])


        #-- have max radius at these bkg locations --#
        def __maxBkgRadius (Bx, By):
            ptRepeatA = np.repeat(np.array([(Bx, By)]), len(otherSrc), axis=0)
            ptRepeatB = np.repeat(np.array([(Bx, By)]), len(cornerLines), axis=0)
            checklist = list(map(__euclideanDist, ptRepeatA, otherSrc))
            checklist = checklist + list(map(__pDistToLine, ptRepeatB, cornerLines))
            checklist.append( __euclideanDist((xC, yC), (Bx, By)) - srcR*(1/3600)/header['CDELT2'] )
            return (min(checklist) - self.gap)  #--leaving a gap around the circle.

        Br1new = __maxBkgRadius(Bx1, By1)
        if Br1new >= Br1:
            Br1 = Br1new
        Br2new = __maxBkgRadius(Bx2, By2)
        if Br2new >= Br2:
            Br2 = Br2new

        #-- check for distance between the two bkg circles --#
        if __euclideanDist((Bx1, By1), (Bx2, By2)) < Br1+Br2:
            if Br1 >= Br2:
                Br1 = __euclideanDist((Bx1, By1), (Bx2, By2)) - Br2 - self.gap
            else:
                Br2 = __euclideanDist((Bx1, By1), (Bx2, By2)) - Br1 - self.gap

        Br = np.array([Br1, Br2])*3600*header['CDELT2']  #--convert to arcsec
    
        #-- plot the randomly selected background points --#
        for i in range(4):
            for Bx, By in zip([Bx1, Bx2], [By1, By2]):
                __pDistToLine((Bx, By), cornerLines[i], ax=axes[3])
                axes[3].plot(Bx, By, 'd', markersize=10, color='darkgreen', markeredgecolor='black', markeredgewidth=0.2)
    
        #-- get 1st background circle --#
        bCirc1 = plt.Circle((Bx1, By1), Br1, alpha=0.50)
        axes[3].add_artist(bCirc1)
    
        #-- get 2nd background circle --#
        bCirc2 = plt.Circle((Bx2, By2), Br2, alpha=0.50)
        axes[3].add_artist(bCirc2)

        #-- zoom to the overlap --#
        pad = 10
        axes[3].set_xlim([min(yiOverlap)-pad, max(yiOverlap)+pad])
        axes[3].set_ylim([min(xiOverlap)-pad, max(xiOverlap)+pad])
        self.axlims['xlim'] = [min(yiOverlap)-pad, max(yiOverlap)+pad]
        self.axlims['ylim'] = [min(xiOverlap)-pad, max(xiOverlap)+pad]
    
        #-- save/show final figure --#
        if self.saveFig:
            plt.savefig(workdir+'/'+'detected_overlap_'+obsID+'.png', dpi=600)
        if self.showFig:
            plt.show()
        plt.close()

        #-- return the background circles --#
        Bx1, By1, Bx2, By2 = np.array([Bx1, By1, Bx2, By2])*header['CDELT2L']  #--coords in Sky coords.
        print('\nFirst Background circle: ', Bx1, By1, Br[0])
        print('Second Background circle: ', Bx2, By2, Br[1])

        self.bCircle1 = [Bx1, By1, Br1*header['CDELT2L']]
        self.bCircle2 = [Bx2, By2, Br2*header['CDELT2L']]
        return True


    ##-- function to create output images --##
    def createOutputImages (self):
        """
        Creates two JPEG output images, one showing the other sources and 
        the other showing the Source and the Background Circles.
        'PNMOS12_image_full' is used as the backdrop.

        :Input: findOverlap object with all previous functions already run.
        :Output: Saved images in the workdir of the obsID.
        """
        workdir, header = self.workdir, self.header
        obsID = self.obsID
        otherSrc = np.array(self.otherSrc)

        if self.isSmallMode:
            img = self.imgPNMOS
        else:
            img = self.overlap

        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        ax.imshow(img, cmap='afmhot', origin='lower')
        ax.set_xlim(self.axlims['xlim'])
        ax.set_ylim(self.axlims['ylim'])

        #-- plot the Source circle --#
        #srcCoords = np.array(self.srcCoords)/header['CDELT2L']  #--convert to log scale pixel.
        srcCoords = self.correctSrc/header['CDELT2L']
        xC, yC = srcCoords
        """
        if self.isSmallMode:
            xC, yC = srcCoords
        else:
            yC, xC = srcCoords   #--the axes are opposite here.
        """
        srcR = self.srcR/header['CDELT2L']  #--converting from log pixel scale.

        ax.plot(xC, yC, '*', markersize=15, color='gray', markeredgecolor='black', markeredgewidth=0.2, label='src')
        ax.add_artist( plt.Circle((xC, yC), srcR, alpha=0.50) )
        ax.add_artist( plt.Circle((xC, yC), srcR, color='blue', fill=False) )

        #-- plot the other sources --#
        ax.plot(otherSrc[:, 0], otherSrc[:, 1], 'p', markersize=10, color='gray', markeredgecolor='black', markeredgewidth=0.2, \
                label='others')
        #for other in otherSrc:
        #    ax.add_artist( plt.Circle((other[0], other[1]), srcR, color='blue', fill=False) )

        #-- plot the background circles --#
        for i, bCircle in enumerate([self.bCircle1, self.bCircle2]):
            Bx, By = bCircle[0], bCircle[1]
            Bx, By = np.array([Bx, By])/header['CDELT2L']
            Br = bCircle[2]/header['CDELT2L']  #--converting from log pixel scale.
            #print(Bx, By, Br, bCircle[2])
            if i==0:
                ax.plot(Bx, By, 'd', markersize=10, color='darkgreen', markeredgecolor='black', markeredgewidth=0.2, label='bkg')
            else:
                ax.plot(Bx, By, 'd', markersize=10, color='darkgreen', markeredgecolor='black', markeredgewidth=0.2)
            ax.add_artist( plt.Circle((Bx, By), Br, alpha=0.50) )
            ax.add_artist( plt.Circle((Bx, By), Br, color='green', fill=False) )

        if self.isSmallMode:
            ax.set_title('small mode', fontsize=12, loc='left')
        ax.set_title(obsID, fontsize=12, loc='right')
        plt.legend(loc='upper right')
        plt.tight_layout()

        #-- save/show final figure --#
        if self.saveFig:
            plt.savefig(workdir+'/'+'source_background_circles.png', dpi=600)
        if self.showFig:
            plt.show()
        plt.close()
        
        return True


    ##-- function to find some background circles --##
    def backgroundCircles_indi (self, axes, inst):
        """
        Automatically finds Background Circles in the instrument data passed.

        :Input: 
            `inst_CCD##_image.fits` files containing the extracted image data of the instrument.
            ``srcCoords`` dictionary containing coords of the main Source.
            ``otherSrc`` dictionary containing coords of other Sources.

        :Requires: Call to `epicObj.find_otherSources_indi` function.

        Args:
            axes (obj): `matplotlib.pyplot.axes` object for plotting the images.
            inst (str): name of the instrument of use.

        :Output: Coordinates and Radii of the Background Circles.

        NOTES
        -----
        ``small-mode`` MOS checking is invalid here since overlap is not been detected.
        """
        workdir, obsID, instCCD = self.workdir, self.obsID, self.pnCCD
        name = inst + '_CCD' + instCCD
        
        #-- open the images --#
        imgInst, header = self.__openImg(workdir+'/'+name+'_image.fits', header=True)

        srcCoords, otherSrc = self.srcCoords, np.array(self.otherSrc)
    
        ##-- distance between two points --##
        def __euclideanDist (P1, P2):
            """
            Calculates the Euclidean distance between two cartesian points.
        
            :Input: Coordinates of the two points in form of tuples.
            :Output: The distance between the points, float.
        
            Note that the Source coordinates are reversed, so that has to be taken 
            into account when passing the values.
            """
            x1, y1 = P1
            x2, y2 = P2
            return np.sqrt( (x2-x1)**2 + (y2-y1)**2 )

        ##-- distance of a point from a line --##
        def __pDistToLine (P0, line, ax=None):
            """
            Calculates the perpendicular distance of point from a line.
        
            :Input: Two tuples P0 and line. P0 contains the x and y coordinates of the point and
                   line contains the coordinates of the start and end points of the line.
            :Output: The perpendicular distance of the point from the line, float.
        
            First the intersection point Pi of the perpendicular from P0 to the line is calculated.
            The output is then the EuclideanDist between P0 and Pi.
        
            Note that the corner points of the overlap are also reversed in coordinates.
            So, the start and end points of the line have reversed x and y.
            """
            x0, y0 = P0
            x1, y1, x2, y2 = line
        
            m = (y2-y1)/(x2-x1)
            mP = -1/m
        
            xI = ((y1-y0) + mP*x0 - m*x1) / (mP-m)
            yI = y0 + mP * (xI-x0)
        
            if ax!=None:
                ax.plot([x1, x2], [y1, y2], '-', color='blue')
                ax.plot(xI, yI, 'p', markersize=10, color='red')
                ax.plot([x0, xI], [y0, yI], '-', color='violet')
            
            return __euclideanDist(P0, (xI, yI))
        
    
        #-- find corners of the overlap --#
        xiOverlap, yiOverlap = self._detect_individual(overlap, wSize=4, percent=0.05, ax=axes[3], mesh=False)
    
        #-- coordinates of all the points in the overlap --#
        xxOverlap, yyOverlap = polygon(yiOverlap, xiOverlap)
        #axes[3].plot(xxOverlap, yyOverlap, '.', color='pink')
       
        #-- lines made by the overlap corners --#
        cornerLines = []
        for i in range(4):
            if i==3:
                cornerLines.append( (yiOverlap[i], xiOverlap[i], yiOverlap[0], xiOverlap[0]) )
            else:
                cornerLines.append( (yiOverlap[i], xiOverlap[i], yiOverlap[i+1], xiOverlap[i+1]) )
    
        #-- the source coordinates --#
        if srcCoords!=None:
            srcCoords = np.array(srcCoords)/header['CDELT2L']
            xC, yC = srcCoords
            """
            if self.isSmallMode:
                xC, yC = srcCoords
            else:
                yC, xC = srcCoords   #--the axes are opposite here.
            """
            print('\nThe Source coordinate value in log scale pixel is: {}, {}'.format(xC, yC))
        else:
            xC, yC = int(imgPNMOS.shape[0]/2), int(imgPNMOS.shape[1]/2)
        #axes[3].plot(xC, yC, '*', markersize=20, color='white', markeredgecolor='black', markeredgewidth=0.2)

        #-- correct the position of main Source --#
        srcRepeat = np.repeat(np.array([(xC, yC)]), len(otherSrc), axis=0)
        oSrcDist = np.array( list(map(__euclideanDist, otherSrc, srcRepeat)) )
        correctSrc = otherSrc[np.where(oSrcDist==min(oSrcDist))]
        xC, yC = correctSrc[0][0], correctSrc[0][1]
        self.correctSrc = np.array([xC, yC])*header['CDELT2L']
        axes[3].plot(xC, yC, '*', markersize=20, color='white', markeredgecolor='black', markeredgewidth=0.2)

        #-- plot other sources --#
        axes[3].plot(otherSrc[:, 0], otherSrc[:, 1], 'p', markersize=10, color='white', markeredgecolor='black', markeredgewidth=0.2)
    
        #-- threshold distance from the Source --#
        dSrcThreshold = 70  #--value in arcsec.
        if self.srcThreshold != None:
            srcThreshold = float(self.srcThreshold)
            if srcThreshold <= dSrcThreshold:
                dSrcThreshold = srcThreshold
            else:
                print(f'\nProvided threshold value exceed {dSrcThreshold} arcsec. Defaulting to this value.')

        DSource = dSrcThreshold*(1/3600)/header['CDELT2']
        print('\nThe distance between random Background circle center and the Source is set at {} arcsec.'.format(dSrcThreshold))

        #-- radius of the background circles --#
        Br = np.array([dSrcThreshold/2, dSrcThreshold/2 - 10])  #--value in arcsec.
        Br1, Br2 = Br*(1/3600)/header['CDELT2']
        if self.bkgR!=None:
            bkgR = np.array(self.bkgR).astype(float)
            if sum(bkgR) <= sum(Br):
                Br = bkgR
                Br1, Br2 = bkgR*(1/3600)/header['CDELT2']
            else:
                print('\nProvided background circle radii are more than half of threshold distance.')
        #print(Br1, Br2)
    
        #-- dSrcThreshold gives the Source circle radius --#
        srcR = dSrcThreshold/2
        srcToLines = [__pDistToLine((xC, yC), cornerLine) for cornerLine in cornerLines]  #--output is in pixels.
        if min(srcToLines) < srcR*(1/3600)/header['CDELT2']:
            print('\nNote: The Source is very near the overlap edge. Source circle area maybe small.')
            srcR = min(srcToLines)*3600*header['CDELT2']  #--convert to arcsec

        if self.srcR != None:
            srcR_given = float(self.srcR)   #--convert from arcsec.
            if srcR_given <= srcR:
                srcR = srcR_given
            else:
                print(f'\nsrcCircRadius value of {np.round(srcR_given, 2)} arcsec provided exceed the distance to the nearest overlap border. Defaulting to this distance.')

        print(f'The Source radius is set at {np.round(srcR, 2)} arcsec.')
        self.srcR = (srcR*(1/3600)/header['CDELT2'] )*header['CDELT2L']  #--convert pixels to log pixels.

        #-- shortlisting random points --#
        xUse, yUse, xUseN, yUseN = [], [], [], []
        minDistToOtherSrc = Br1
        Bx1, By1 = None, None
        for ptX, ptY in zip(xxOverlap, yyOverlap):
            ptRepeat0 = np.repeat(np.array([(ptX, ptY)]), len(cornerLines), axis=0)
            pDistList = list(map(__pDistToLine, ptRepeat0, cornerLines))
            pDistList.sort()

            #-- check distance from main Source and cornerLines --#
            if __euclideanDist((ptX, ptY), (xC, yC)) > (DSource + self.gap) and pDistList[0] > (Br1 + self.gap):
                    xUse.append(ptX)
                    yUse.append(ptY)

                    #-- check distance from other sources --#
                    ptRepeat = np.repeat(np.array([(ptX, ptY)]), len(otherSrc), axis=0)
                    distList = list(map(__euclideanDist, ptRepeat, otherSrc))
                    distList.sort()
                    if distList[0] > DSource:
                        xUseN.append(ptX)
                        yUseN.append(ptY)
                        if distList[0] >= (minDistToOtherSrc + self.gap):
                            Bx1, By1 = ptX, ptY
 
        #print(len(xUseN), len(yUseN))
        
        ##-- find random background points --##
        def __backgroundPt (xToUse=xUseN, yToUse=yUseN, bkgPt=None, ax=None):
            """ 
            Internal function to select random constrained points.
            At first, the arrays (xUseN, yUseN) are used to find the random points.

            If the number of points in these is small and two random points atleast the 
            threshold distance away from eachother are not found, then the arrays (xUse, yUse)
            are used. 
            """
           
            message = f'\nNOTE:There are too many other sources for obsID {obsID}, such that at the threshold given, no background circles are possible.' \
                        +'\nKindly redo the product extraction for this one later with some changed parameters.' \
                        +'\nFor now the other sources are not taken into account for obtaining the background circles.\n'
            
            #-- check if the arrays are empty --#
            if len(xToUse)<2:
                print(message)
                return __backgroundPt(xToUse=xUse, yToUse=yUse, ax=ax)

            #-- first random point --#
            if bkgPt is None:
                idx1 = np.random.choice( range(len(xToUse)) )
                Bx1, By1 = xToUse[idx1], yToUse[idx1]
            else:
                Bx1, By1 = bkgPt

            #-- second background point --#
            totalBr = Br1 + Br2
            Bx2, By2 = None, None
            for ptX, ptY in zip(xToUse, yToUse):
                dist = __euclideanDist((Bx1, By1), (ptX, ptY))
                if dist >= totalBr:
                    totalBr = dist
                    Bx2, By2 = ptX, ptY

            if Bx2 is not None:
                axes[3].plot(xToUse, yToUse, '.', color='cyan', alpha=0.25)
                return [Bx1, By1, Bx2, By2]
                
            #-- when other sources are too many --#
            else:
                print(message)
                return __backgroundPt(xToUse=xUse, yToUse=yUse, ax=ax)

        if Bx1 is not None:
            Bx1, By1, Bx2, By2 = __backgroundPt(bkgPt=(Bx1, By1), ax=axes[3])
        else:
            Bx1, By1, Bx2, By2 = __backgroundPt(ax=axes[3])


        #-- have max radius at these bkg locations --#
        def __maxBkgRadius (Bx, By):
            ptRepeatA = np.repeat(np.array([(Bx, By)]), len(otherSrc), axis=0)
            ptRepeatB = np.repeat(np.array([(Bx, By)]), len(cornerLines), axis=0)
            checklist = list(map(__euclideanDist, ptRepeatA, otherSrc))
            checklist = checklist + list(map(__pDistToLine, ptRepeatB, cornerLines))
            checklist.append( __euclideanDist((xC, yC), (Bx, By)) - srcR*(1/3600)/header['CDELT2'] )
            return (min(checklist) - self.gap)  #--leaving a gap around the circle.

        Br1new = __maxBkgRadius(Bx1, By1)
        if Br1new >= Br1:
            Br1 = Br1new
        Br2new = __maxBkgRadius(Bx2, By2)
        if Br2new >= Br2:
            Br2 = Br2new

        #-- check for distance between the two bkg circles --#
        if __euclideanDist((Bx1, By1), (Bx2, By2)) < Br1+Br2:
            if Br1 >= Br2:
                Br1 = __euclideanDist((Bx1, By1), (Bx2, By2)) - Br2 - self.gap
            else:
                Br2 = __euclideanDist((Bx1, By1), (Bx2, By2)) - Br1 - self.gap

        Br = np.array([Br1, Br2])*3600*header['CDELT2']  #--convert to arcsec
    
        #-- plot the randomly selected background points --#
        for i in range(4):
            for Bx, By in zip([Bx1, Bx2], [By1, By2]):
                __pDistToLine((Bx, By), cornerLines[i], ax=axes[3])
                axes[3].plot(Bx, By, 'd', markersize=10, color='darkgreen', markeredgecolor='black', markeredgewidth=0.2)
    
        #-- get 1st background circle --#
        bCirc1 = plt.Circle((Bx1, By1), Br1, alpha=0.50)
        axes[3].add_artist(bCirc1)
    
        #-- get 2nd background circle --#
        bCirc2 = plt.Circle((Bx2, By2), Br2, alpha=0.50)
        axes[3].add_artist(bCirc2)

        #-- zoom to the overlap --#
        pad = 10
        axes[3].set_xlim([min(yiOverlap)-pad, max(yiOverlap)+pad])
        axes[3].set_ylim([min(xiOverlap)-pad, max(xiOverlap)+pad])
        self.axlims['xlim'] = [min(yiOverlap)-pad, max(yiOverlap)+pad]
        self.axlims['ylim'] = [min(xiOverlap)-pad, max(xiOverlap)+pad]
    
        #-- save/show final figure --#
        if self.saveFig:
            plt.savefig(workdir+'/'+'detected_overlap_'+obsID+'.png', dpi=600)
        if self.showFig:
            plt.show()
        plt.close()

        #-- return the background circles --#
        Bx1, By1, Bx2, By2 = np.array([Bx1, By1, Bx2, By2])*header['CDELT2L']  #--coords in Sky coords.
        print('\nFirst Background circle: ', Bx1, By1, Br[0])
        print('Second Background circle: ', Bx2, By2, Br[1])

        self.bCircle1 = [Bx1, By1, Br1*header['CDELT2L']]
        self.bCircle2 = [Bx2, By2, Br2*header['CDELT2L']]
        return True


    ##-- main function --##
    def main (self, inst=None):
        """
        The main function to sequentially run the functions in `findOverlap` class.

        Returns:
            bCircle1 (tuple): (x, y, r) coords and radius of first background circle.
            bCircle2 (tuple): (x, y, r) coords and radius of second background circle.

            correctSrc (tuple): corrected (x, y) coords of the main Source.
                The original main Source coords found by the `epicObj.findSourceCCD` 
                function using the XMMSAS ``ecoordconv`` utility are a little bit shifted
                from what the main Source likely is as found from ``edetect_chain`` while
                finding other Sources. So, the main Source coords are changed to the coords
                of the nearest Source to it in the ``otherSrc`` dictionary.

            srcR (float): radius of the Source circle or outer radius of the Source Annulus
                for the Piled-up cases, although a full circle is detected here.

            isSmallMode (bool): True, if absolute difference between flattened MOS12 and 
                PNMOS12 image arrays is less than 500. False otherwise.
                See issue #28 in the xmmPipeline repository in this regard, 
                https://github.com/Suyog7130/xmmPipeline/issues/28

        NOTES
        -----
        All returned coordinates are in **Physical Units**.
        These (x, y) are termed as (**RAWX**, **RAWY**) in DS9.
        The radius values are in **Sky Coordinates**, which I think is the name XMMSAS gives
        to **Physical Units**. Each Sky Coord unit equals **0.05 arcsec**.
        So, a radius of 640 equals 32 arcsec.
        See issue #24 for more details, https://github.com/Suyog7130/xmmPipeline/issues/24
        """
        fig, axes = plt.subplots(1, 4, figsize=(20, 5))

        #-- if finding srcBkg circs in overlap --#
        if inst is None:
            self.obtainOverlap(axes=axes)
            self.backgroundCircles(axes=axes)
            self.createOutputImages()
            return (self.bCircle1, self.bCircle2, self.correctSrc, self.srcR, self.isSmallMode)
        
        #-- for individual instrument --#
        self.backgroundCircles(axes=axes)
        self.createOutputImages()
        return (self.bCircle1, self.bCircle2, self.correctSrc, self.srcR)



if __name__=="__main__":

    maindir = '/media/suyog/DATA/xmm_obs'
    obsIDs = ['0770980601', '0810200701', '0831790201']
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--obsIDs', nargs='+', action='store', default=obsIDs)
    
    args = parser.parse_args()
    obsIDs = args.obsIDs
    
    start = time.time()
    for obsID in obsIDs:
        print('\nWorking for obsID {}.'.format(obsID))
        workdir = maindir+'/'+str(obsID)+'/work'

        findOverlapObj = findOverlap(workdir=workdir, obsID=obsID, saveFig=False)
        bCircle1, bCircle2, srcR, isSmallMode = findOverlapObj.main()
        print(bCircle1, bCircle2, srcR, isSmallMode, sep='\n')
    end = time.time()
    print(f'Time taken for execution: {np.round(end-start, 5)} sec')
    
    
    
    
    
    
#################### End of Program #########################
#############################################################
    
    
    
    
    
    
    
    
    
    
    
    
    

