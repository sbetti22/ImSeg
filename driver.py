## written by Sarah Betti 2019
## updated 5 November 2019

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.wcs import WCS

import glob
import os
import time
import re

from getrot import *
from coreSelection import *
from getBlobs import *
from noiseblobs import *
plt.ion()

############## CHANGEABLE PARAMETERS ##################

#PATH TO SIGNAL/WEIGHT/SN/PSF MAPS
#signal  = 'coadded_signal.fits'
#weight  = 'coadded_weight.fits'
#sn      = 'coadded_sn.fits'
#psf     = 'coadded_psf.fits'
signal  = 
weight  = 
sn      = 
psf     = 

#PATH TO NOISE REALIZATION FILES.  DO NOT INCLUDE NAMES OF FILES!
#noisepath = '/path/to/noise/realization/' OR None
noisepath = 
 
#OBSERVATION/TARGET INFORMATION
#wavelength = wavelength of observations in microns
#pcdist = distance to target in parsecs       
#sn_thresh = S/N threshold for image segmenataion
#sim = define whether or not maps are simulations (True) or observations (False)
wavelength = 
pcdist = 
sn_thresh = 
sim = 

#PATH TO TEMPERATURE/COLUMN DENSITY/YSO LISTS
#path to temperature maps
#temp = 'temperature.fits'
temp = 

#path to column density map
#coldens = 'coldens.fits'
coldens = 

#load YSO list
#ysoradec = np.loadtxt('YSOs.txt')
ysoradec = 

#LOCATION TO SAVE OUTPUT FILES
#outbase='/path/to/output/files' 
outbase = 

#PARTS OF REDUCTION TO RUN
run_getblobs = True
run_getnoiseblobs = True 
run_coreselection = True


############### END CHANGEABLE PARAMETERS ########################


############### RUN IMSEG ########################
#make output file folder
if not os.path.exists(outbase):
    os.makedirs(outbase)

t1 = time.time()

print('starting core extraction')
print(f'path to signal files: {signal}')
print(f'path to noise files: {noisepath}')  
print(f'path to output files: {outbase} \n')
print('---------------o---------------\n')

mb = getBlobs('cores',pcdist, signal,weight,sn,psf,outbase, temp=temp, temp_hdr=temp, sim=sim)

if run_getblobs: 
    print('starting getBlobs\n')
    mb.nsigma = 5 #5
    mb.nmxsigma = 3.5 #3.5
    mb.areacut = -100
    mb.covlim = 0.4
    mb.contrast = 0.001
    #set number of connected pixels based on distance to target
    npix = int((206265*0.04)/pcdist)
    mb.npixels = npix 

    mb.getBlobCatalogue(sn_thresh,ysoradec, plot_blobs=True)
    print('---------------o---------------\n')
    
if run_getnoiseblobs:
    print('starting getNoiseBlobs\n')
    runNoiseBlobs(pcdist, noisepath, signal, psf, outbase, sn_thresh,sim, cleanup=False)
    print('---------------o---------------\n')
    
fwhm_beam = mb.FWHM()
pickle_file = mb.bloblist_name()
noise_pickle_file = noisebloblist_name()
cs = coreSelection(outbase, wavelength, pcdist, coldens, pickle_file, noise_pickle_file=noise_pickle_file)
if run_coreselection:
    print('starting runCoreSelection\n')
    good_indices, goodmap = cs.runCoreSelection()
    cs.saveGoodCores(fwhm_beam, good_indices=good_indices)
    cs.plot_ratio('sim', sim=sim, good_indices=good_indices, goodmap=goodmap)
    print('---------------o---------------\n')
    
t2 = time.time()
print('finished total time: ', (t2-t1)/60., ' min')
