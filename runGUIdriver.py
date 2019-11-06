## written by Sarah Betti 2019
## updated 5 November 2019

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.wcs import WCS
from sys import argv

import glob
import os
import time

from getrot import *
from coreSelection import *
from getBlobs import *
from noiseblobs import *

plt.ion()

for i in np.arange(len(argv)):
    if argv[i] == 'None':
        argv[i] = None
    if (argv[i] == 'False') or (argv[i] == '0'):
        argv[i] = False
    if (argv[i] == 'True') or (argv[i] == '1'):
        argv[i] = True
    

script, name, signal, sn, weight, psf, noisepath, wavelength, pcdist, sn_thresh, sim, temp, wcs, coldens, yso, outbase, run_getblobs, run_getnoiseblobs, run_coreselection, ysora, ysodec, instrument = argv

print('---------')

sn_thresh = float(sn_thresh)
wavelength = int(wavelength)
pcdist = int(pcdist)
ysora = int(ysora)
ysodec = int(ysodec)


if not os.path.exists(outbase):
    os.makedirs(outbase)

yso_list = np.loadtxt(yso, comments='#')
yso_ra = yso_list[:,ysora]
yso_dec = yso_list[:,ysodec]
ysoradec = np.c_[yso_ra, yso_dec]

t1 = time.time()

print(f'starting core extraction for {name}')

if signal:
    print(f'path to signal files: {signal}')
if noisepath is not None:
    print(f'path to noise files: {noisepath}')

print(f'path to output files: {outbase} \n')

print('---------------o---------------\n')

mb = getBlobs(name,pcdist, signal,weight,sn,psf,outbase, temp=temp, temp_hdr=wcs, sim=sim)
if run_getblobs: 
    print('starting getBlobs\n')
    mb.nsigma = 5 #5
    mb.nmxsigma = 3.5 #3.5
    mb.areacut = -100
    mb.covlim = 0.4
    npix = int((206265*0.04)/pcdist)
    mb.npixels = npix
    mb.contrast = 0.001

    mb.getBlobCatalogue(sn_thresh,ysoradec, plot_blobs=True)
    print('---------------o---------------\n')

if run_getnoiseblobs:
    print('starting getNoiseBlobs\n')
    runNoiseBlobs(pcdist, noisepath, signal, psf, outbase, sn_thresh,sim, cleanup=False)
    print('---------------o---------------\n')

fwhm_beam = mb.FWHM()
pickle_file = mb.bloblist_name()

noise_pickle_file = noisebloblist_name()
if not os.path.exists(f'{outbase}/{noise_pickle_file}'):
    noise_pickle_file = None
cs = coreSelection(outbase, wavelength, pcdist, coldens, pickle_file, noise_pickle_file=noise_pickle_file)

if run_coreselection:
    print('starting runCoreSelection\n')
    good_indices, goodmap = cs.runCoreSelection()
    
    cs.saveGoodCores(fwhm_beam, good_indices=good_indices)
    cs.plot_CMF()
    cs.plot_ratio(instrument, sim=sim, good_indices=good_indices, goodmap=goodmap)
    
    print('---------------o---------------\n')

t2 = time.time()
print('finished total time: ', (t2-t1)/60., ' min')
