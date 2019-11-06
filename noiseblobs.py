## written by Sarah Betti 2019
## updated 5 November 2019

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.wcs import WCS

import time
import os
import glob

from getBlobs import *


def setup_NoiseBlobs(ncfile, signalheader, psfdata):
    #go through noise.nc files and create fits files
    f1 = ncfile
    
    if (os.path.isfile(f1.split('.nc')[0] + '_weight.fits') == False) or (os.path.isfile(f1.split('.nc')[0] + '_signal.fits') == False):
        print(f'creating signal, weight, sn noise maps for noise map: {f1.split("/")[-1]}')
        os.system('fitswriter ' + f1)
        
    #create sn and psf files
    if os.path.isfile(f1.split('.nc')[0] + '_sn.fits') == False:  
        fil = f1.split('.nc')[0]
        imwt = fits.getdata(fil + '_weight.fits')
        imsig = fits.getdata(fil + '_signal.fits')
        imsn = imsig*np.sqrt(imwt)
        
        #determine if any weights are less than 1.  if so, make the S/N at that location 0 and not negative. 
        w0 = np.where(imwt < 1.)
        n0 = len(w0[0])
        if n0 > 0:
            imsn[w0] = 0.
            
        fits.writeto(fil + '_sn.fits', imsn, signalheader, overwrite=True)
        
    if os.path.isfile(f1.split('.nc')[0] + '_psf.fits') == False:
        fil = f1.split('.nc')[0]
        fits.writeto(fil + '_psf.fits', psfdata, signalheader, overwrite=True)
        

def noiseMaps(noisepath, psf, signal):
    '''
    find all noise realization maps  
    
    Parameters
    ----------
    path : str
        location of coadd.nc
        
    Returns
    ---------
    globbed list of path location of all noise realization maps within /noise_maps/ folder.  
    
    '''

    f = glob.glob(noisepath + '/noise*.nc')  
    
    if len(f) > 0:
        ## reorder files so that they are in a natural/human ordering
        noise_number = [int(re.findall('\d+', i.split('/')[-1])[0]) for i in f] #pull out integer number of noise map 
        noise_num_ind = np.argsort(noise_number)  #reorder based on noisemap value
        f = np.array(f)[noise_num_ind] #sort initial list based on reordering
        
        # get coadded signal and psf
    psfdata = fits.getdata(psf)
    signalheader = fits.getheader(signal)
    
    for ncfile in f:
        setup_NoiseBlobs(ncfile, signalheader, psfdata)
        
    return f


def blobListSav(outbase):
    '''
    find all noise bloblist sav file
        
    Returns
    ---------
    globbed list of path location of all noise realization bloblist save files within /noise_maps/ folder.  
    '''
    f = glob.glob(outbase + '/bloblist_noise*.pkl')
    
    ## reorder files so that they are in a natural/human ordering
    noise_number = [int(re.findall('\d+', i.split('/')[-1])[0]) for i in f] #pull out integer number of noise map 
    noise_num_ind = np.argsort(noise_number)  #reorder based on noisemap value
    f = np.array(f)[noise_num_ind] #sort initial list based on reordering
    
    return f


def noise_maps(noisepath):
    signal = noisepath + '_signal.fits'
    weight = noisepath+ '_weight.fits'
    sn = noisepath + '_sn.fits'
    psf = noisepath + '_psf.fits'
    
    return signal, weight, sn, psf

        
def runNoiseBlobs(pcdist, noisepath, signal, psf, outbase, sn_thresh, sim, cleanup=False):
    '''
    determine noise blobs in each noise map by running each noise .nc map through self.runMacanaBlobs() and save resulting properties in pickled pandas dataframe 
    
    Parameters
    ----------
    coaddpath : str
        location of coadd.nc
    noisepath : str
        location of noise*.nc files
    pcdist : float
        distance to source in pc
    cleanup : boolean 
        if True - delete all .fits files for noise*.nc 
        default is False
    
            
    Returns
    ---------
    saves pickled pandas dataframe of all noise core properties
    '''
    
    #get noise*.nc fies
    

    f1 = noiseMaps(noisepath, signal, psf)
    
    #run macana_blobs on noise images
    for i in np.arange(len(f1)):
        fil = f1[i].split('.nc')[0]
        num = fil.split('/noise')[-1]
        
        noisesignal, noiseweight, noisesn, noisepsf = noise_maps(fil)
        
        if os.path.isfile(outbase + '/bloblist_noise' + str(num) + '.pkl') == False:
            
            gb = getBlobs('noise'+num, pcdist, noisesignal,noiseweight,noisesn,noisepsf,outbase, sim=sim)
            
            npix = int((206265*0.04)/pcdist)
            gb.npixels = npix 
            gb.getBlobCatalogue(sn_thresh, ysoradec = None, plot_blobs=False)

            print()   
        else:
            print('noise map ' + str(num) + ' have already been segmented: output files in: ' + outbase + '/bloblist_noise' + str(num) + '.pkl')

        if cleanup == True:
            os.system('rm ' + fil + '_*.fits')

    #open pickle files from macana blobs and combine into one noise blobs file
    f = blobListSav(outbase)
    print()
    print(f'number of bloblist.pkl noise files: {len(np.array(f))}')
    
    total_noise = []
    for i in f:
        c = pd.read_pickle(i)
        total_noise.append(c)
    total_noise = pd.concat(total_noise, axis=0)
    
    df_total_noise = pd.DataFrame(data = total_noise)
    df_total_noise = df_total_noise.reset_index()
    df_total_noise['numnoise'] = np.ones_like(df_total_noise['ra'])*len(f)
    
    df_total_noise.to_pickle(outbase + '/bloblist_total_noise.pkl')
    
    print(f'saved final noisebloblist to: {outbase}/bloblist_total_noise.pkl')
    
    
def noisebloblist_name():
    return 'bloblist_total_noise.pkl'


    
    


