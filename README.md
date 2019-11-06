
# Introduction

ImSeg is an molecular cloud core detection and characterization utility which uses the phoutils package `photutils.segmentation` to detect and deblend point-like and extended sources in both observational data and simulations.  

It is written specifically for LMT AzTEC and TolTEC data analysis, but can be applied to other data.  

After detecting core candidates, ImSeg then measures the characteristics of each core, including ra/dec, total and peak flux, area, fwhm major/minor axes, and mass.  All core candidates are then put through a "goodness" test to determine the final core catalogue and core mass function (CMF).


# Packages required
The following packages must be installed for ImSeg to work:
- `numpy`
- `pandas`
- `astropy`
- `glob`
- `os`
- `time`
- `re`
- `matplotlib`
- `photutils`
- `scipy`





Before running ImSeg, several changes need to be made to the `photutils.segmentation` function `deblend_sources()` in deblend.py.  These changes are necessary so that the multi-thresholding is set to 1 sigma levels which allows for blobs of different S/N to be segmented evenly.

Within deblend.py: 

in `deblend_sources()` function make the following changes:

change: 

> ```def deblend_sources(data, segment_img, npixels, filter_kernel=None, labels=None, nlevels=32, contrast=0.001, mode='exponential', connectivity=8, relabel=True):```

to: 

> ```def deblend_sources(data, data_byte,segment_img, npixels, filter_kernel=None,labels=None, contrast=0.001, mode='exponential', connectivity=8, relabel=True):```

--------------------

change: 

>```for label in labels:
    source_slice = segment_img.slices[label - 1]
    source_data = data[source_slice]
    source_segm = SegmentationImage(np.copy(segment_img.data[source_slice]))
    source_segm.keep_labels(label)    # include only one label
    source_deblended = _deblend_source(source_data, source_segm, npixels, nlevels=nlevels, contrast=contrast, mode=mode, connectivity=connectivity) ```




to: 

> ```for label in labels:
    source_slice = segment_img.slices[label - 1]
    source_data = data[source_slice]
    source_data_byte = data_byte[source_slice]
    nlevels = 255 - np.unique(source_data_byte)
    nlevels = nlevels[::-1]
    source_segm = SegmentationImage(np.copy(segment_img.data[source_slice]))
    source_segm.keep_labels(label)    # include only one label
    source_deblended = _deblend_source(source_data, source_data_byte, source_segm, npixels, nlevels=nlevels, contrast=contrast, mode=mode, connectivity=connectivity) ```

--------------------

In ```_deblend_source()``` function, make the following changes:

change: 

>```def _deblend_source(data, segment_img, npixels, nlevels=32, contrast=0.001, mode='exponential', connectivity=8):```
                    
to:                    

>```def _deblend_source(data, data_byte, segment_img, npixels, nlevels, contrast=0.001, mode='exponential', connectivity=8):```

--------------------

change: 
>```if nlevels < 1:
        raise ValueError('nlevels must be >= 1, got "{0}"'.format(nlevels))```


to: 

>``` '''if nlevels < 1:
        raise ValueError('nlevels must be >= 1, got "{0}"'.format(nlevels))''' ```

--------------------

change: 

>```steps = np.arange(1., nlevels + 1)
if mode == 'exponential':
    if source_min == 0:
        source_min = source_max * 0.01
    thresholds = source_min * ((source_max / source_min) **
                                   (steps / (nlevels + 1)))
elif mode == 'linear':
    thresholds = source_min + ((source_max - source_min) /
                                   (nlevels + 1)) * steps
else:
    raise ValueError('"{0}" is an invalid mode; mode must be '
                         '"exponential" or "linear"')```
                         
                         
        
        
to: 

>```if mode=='exponential':
    if source_min == 0:
        source_min = source_max * 0.01
    thresholds = source_min * ((source_max / source_min) ** \
				   (nlevels / (np.max(data_byte) - np.min(data_byte))))
elif mode == 'linear':
    thresholds = source_min + ((source_max - source_min) / \
				   (np.max(data_byte)-np.min(data_byte)))*nlevels
else:
    raise ValueError('"{0}" is an invalid mode; mode must be '
                         '"exponential" or "linear"')```
                         
                         
--------------------




# Inputs required

To run ImSeg, the following files are needed:
- signal .fits map
- weight .fits map
- signal to noise .fits map
- psf .fits map
- N(H$_2$) .fits map
- (optional) noise realization (.nc or .fits) maps (if .fits, must have signal/weight/s2n/psf maps)
- (optional) temperature .fits map
- (optional) yso ra/dec list



# Setup & Run ImSeg

To run ImSeg, two methods are available.  1) a python GUI (GUI_driver.py) or 2) python script (driver.py)

If using the GUI: input all parameters using the buttons.

If using the script: input all parameters within the CHANGEABLE PARAMETERS comments

The following takes you through setting up the parameters within the python script.

## Set up paths and file locations

###  a. set up paths to the data and no realization files 


```python
signal = '/path/to/signal/map/signal.fits'
weight = '/path/to/weight/map/weight.fits'
sn= '/path/to/sn/map/sn.fits'
psf = '/path/to/psf/map/psf.fits'

#(optional: only path to files.  Include a / at the end)
noisepath = '/path/to/noise/realizations/'
```

### b. set up observation and target information:
- wavelength (in microns)
- distance to object (in pc)
- S/N threshold for core candidacy
- sim: (boolean) whether of not the data is from a simulation


```python
wavelength = 1100 #um
pcdist = 830 #pc
sn_thresh = 2.5
sim = False #(True)
```

### c. set up paths to column density, and optional temperature map


```python
#(optional)
temp = '/path/to/temperature/map/temp.fits'

colden = '/path/to/column/density/map/colden.fits'

```

### d. set up YSO ra/dec list

Two options:
1. ysoradec = 2d numpy.ndarray of the form:
    >`ysoradec = array([[ra1, dec1], [ra2, dec2], ... [ran, decn]])`
2. `ysoradec = None`


```python
ysoradec = np.array([ra_1, dec_1], [ra_2, dec_2])
# OR 
ysoradec = np.loadtxt('/path/to/yso/list.txt')
```

### e. set up the location to save all file produced from Imseg.  Do not include a `/` at the end.


```python
outbase = '/path/to/save/all/output/files'
```

### e. determine which parts of reduction to run:

>`run_getblobs` : runs blob identification and characterization from maps

>`run_getnoiseblobs` : run blob identification and characterization from noise realization maps

>`run_coreselection` : run a core selection process to determine core candidates from all blobs found from getBlobs



```python
run_getblobs = True #(False)
run_getnoiseblobs = True #(False)
run_coreselection = True #(False)
```

These are all the initial changeable parameters necessary to run ImSeg.

## To run ImSeg, in a terminal: 

`$ python driver.py`

or select `RUN` if using the GUI


There are three parts of the reduction process.  We will describe how each one works below. 
## Run `getBlobs`

Now, let's start by detecting all blobs and determining characteristics from .fits data maps provided by the `getBlobs` module: 

First, create an instance of `getBlobs` for the .fits data maps


```python
from getBlobs import *
mb = getBlobs('cores', pcdist, signal, weight, sn, psf, outbase, temp=temp, temp_hdr=temp, sim=sim)
```

The source segmentation and characterization is performed by the `getBlobs.getBlobCatalogue` attribute which uses the `clumpFind` class.  This class calls utilizes the `photutils` image segmentatation functions `detect_sources()` and `deblend_sources()` which use a mixture of mulit-thresholding and watershed to detect blobs.  For more information about `photutils.segmentation` see: https://photutils.readthedocs.io/en/stable/segmentation.html#module-photutils.segmentation 

There are several class variables which can be initialized after creating an instance of `getBlobs`.  Changing these class variables is optional.  If no manual inputs are given, the default values will be used.


```python
mb.nsigma = 5
mb.nmxsigma=3.5
mb.areacut=-100
mb.covlim=0.4
mb.npixels=10
mb.contrast=0.001
```

Let's run `getBlobs.getBlobCatalogue` to find all blobs and their properties.  Several outputs will be saved into the outbase directory created above.  


```python
mb.getBlobCatalogue(sn_thresh, ysoradec, plot_blobs=True)
```

## Run `NoiseBlobs`

Now, let's run the function `runNoiseBlobs()` to detect all blobs from the noise realization maps. If no noise realization maps exist, skip this step.  

`runNoiseBlobs()` will create a new instance of `getBlobs` for each noise realization map to determine all blobs.  All outputs will be saved to the outbase directory.


```python
from noiseblobs import *
runNoiseBlobs(pcdist, noisepath, signal, psf, outbase, sn_thresh, sim, cleanup=False)
```

## Run `coreSelection`

The core candidates are selected using the coreSelection module.  The candidates are selected using a "goodness" test based on the whether the core passes 3 tests:
- S/N threshold
- "good" score minimum threshold
- minimum column density ratio threshold


First, initalize an instance of `coreSelection` for the blobs.  To do this, you need to either input the following values or use your `getBlobs` instance to determine the values:  
- fwhm_beam : fwhm of the beam of the data image
- pickle_file : name of the .pkl file created from `getBlobCatalogue`
- noise_pickle_file (optional) : name of the noise .pkl file created from `runNoiseBlobs()` function


```python
fwhm_beam = mb.FWHM()
pickle_file = mb.bloblist_name()
noise_pickle_file = noisebloblist_name()
```


```python
from coreSelection import *
cs = coreSelection(outbase, wavelength, pcdist, coldens, pickle_file, noise_pickle_file = noise_pickle_file)
```

To get the indices from the pickle_file .pkl table of core candidates and a 2D histogram of the goodness selection, run `coreSelection.runCoreSelection()` 


```python
good_indices, goodmap = cs.runCoreSelection()
```

To make and save a new catalogue of the core candidates (as a .csv file), run `coreSelection.saveGoodCores()`


```python
cs.saveGoodCores(fwhm_beam, good_indices=good_indices)
```

To plot the ratio of N(H$_2$) of two wavelengths (colden vs signal maps) to S/N, run `coreSelection.plot_ratio()`


```python
cs.plot_ratio('Observation', sim=sim, good_indices=good_indices, goodmap=goodmap)
```

To plot the CMF of core candidates, run `coreSelection.plot_CMF()` 


```python
cs.plot_CMF()
```

# Outputs

Running `getBlobs` will produce the following outputs:
- bloblist_core.pkl : pickle file of blob properties
- clumpfind_cores.pdf : pdf of original signal map (left) and blobs found (right)
- core_blobprints.fits : fits file of blobs where each blob is labeled by a different integer value.  Zero is always reserved for the background.

Running `runNoiseBlobs()` will produce the save outputs as `getBlobs` for each noise realization along with bloblist_total_noise.pkl, a pickle file of the combined blob properties from all noise realizations.

Running `coreSelection` will produce the following outputs: 
- core_mass_select.npz : core mass selection criteria
- final_core_selection.csv : catalogue of core properties (if `coreSelection.saveGoodCore()` is run) 
- SNcolumndensity.pdf : plot of ratio of column density vs S/N (if `coreSelection.plot_ratio()` is run)
- goodmap.npz : (x,y,z) to create contours on SNcolumndensity.pdf (if `coreSelection.plot_ratio()` is run)
- SNvsNH2ratio.npz : (x,y) data to create all blobs on SNcolumndensity.pdf (if `coreSelection.plot_ratio()` is run)
- SNvsNH2ratio_good.npz : (x,y) data to create all core candidates on SNcolumndensity.pdf (if `coreSelection.plot_ratio()` is run)
- histogram_cmf.pdf : CMF plot (if `coreSelection.plot_CMF()` is run)

## Read output files

All pickle files can be read with `pandas`: 

    pd.read_pickle('/path/to/pickle/file/file.pkl')

All csv files can be read with `pandas`:

    pd.read_csv('/path/to/csv/file/file.csv')

All numpy save files can be read with `numpy`:

    np.load('/path/to/numpy/savefile/file.npz')
