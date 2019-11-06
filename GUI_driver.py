## written by Sarah Betti 2019
## updated 5 November 2019

import numpy as np
import os

import tkinter as tk
from tkinter import *
from tkinter import ttk, StringVar, IntVar, BooleanVar
from tkinter.filedialog import askopenfilename, askdirectory

from PIL import ImageTk, Image


class Window:       
    def __init__(self,window):  
        
        self.name = StringVar()
        
        self.instrument = StringVar()
        
        self.signal = StringVar()
        self.signalpath = ''
        
        self.weight = StringVar()
        self.weightpath = ''
        
        self.sn = StringVar()
        self.snpath = ''
        
        self.psf = StringVar()
        self.psfpath = ''
        
        self.noisepath = StringVar()
        
        self.temp = StringVar()
        self.temppath = ''
        
        self.wcs = self.temp
        
        self.nh2 = StringVar()
        self.nh2path = ''
        
        self.yso = StringVar()
        self.ysopath = ''
        
        self.racol = IntVar()
        self.deccol = IntVar()
        
        self.outbase = StringVar()

        
        self.wave = IntVar()
        self.pcdist = IntVar()
        self.sn_thresh = IntVar()
        self.sim = BooleanVar()
        
        window.title('ImSeg Reduction')
        
        ttk.Style().configure('red.TButton', background='red')
        
        #name 
        y = 0
        ttk.Label(window, text='Name').grid(row=y)
        ttk.Entry(window, textvariable=self.name).grid(row=y,column=1, ipadx=3, ipady=5, sticky='W')    
        
        y+=1
        ttk.Label(window, text='Instrument').grid(row=y)
        ttk.Entry(window, textvariable=self.instrument).grid(row=y,column=1, ipadx=3, ipady=5, sticky='W')   
        
        y+=1
        ttk.Button(window, text='Get signal\nmap', command = lambda : self.load_signal()).grid(row=y, ipadx=5, ipady=15, sticky="ew")
        
        ttk.Entry(window, textvariable=self.signal, width=70).grid(row=y,column=1, ipadx=1, ipady=1,sticky='W')
        
        ttk.Button(window, text='Open fits', command = lambda : self.process_fits(self.signalpath)).grid(row=y, column=2, ipadx=5, ipady=15, sticky="ew")
        
        
        
        y+=1
        #sn
        ttk.Button(window, text='Get S/N\nmap', command = lambda : self.load_sn()).grid(row=y, ipadx=5, ipady=15, sticky="ew")
        
        ttk.Entry(window, textvariable=self.sn, width=70).grid(row=y,column=1, ipadx=1, ipady=1,sticky='W')
        
        ttk.Button(window, text='Open fits', command = lambda : self.process_fits(self.snpath)).grid(row=y, column=2, ipadx=5, ipady=15, sticky="ew")
        
        #weight
        y+=1
        ttk.Button(window, text='Get weight\nmap', command = lambda : self.load_weight()).grid(row=y, ipadx=5, ipady=15, sticky="ew")
        
        ttk.Entry(window, textvariable=self.weight, width=70).grid(row=y,column=1, ipadx=1, ipady=1,sticky='W')
        
        ttk.Button(window, text='Open fits', command = lambda : self.process_fits(self.weightpath)).grid(row=y, column=2, ipadx=5, ipady=15, sticky="ew")
        
        #psf
        y+=1
        ttk.Button(window, text='Get psf\nmap', command = lambda : self.load_psf()).grid(row=y, ipadx=5, ipady=15, sticky="ew")
        
        ttk.Entry(window, textvariable=self.psf, width=70).grid(row=y,column=1, ipadx=1, ipady=1,sticky='W')
        
        ttk.Button(window, text='Open fits', command = lambda : self.process_fits(self.psfpath)).grid(row=y, column=2, ipadx=5, ipady=15, sticky="ew")
        
        
        #noise path
        y+=1
        ttk.Button(window, text='Get noisepath', command = lambda : self.load_noise_path()).grid(row=y, ipadx=5, ipady=15, sticky="ew")
        
        ttk.Entry(window, textvariable=self.noisepath, width=70).grid(row=y,column=1, ipadx=1, ipady=1,sticky='W')
        
        y+=1
        ttk.Label(window, text='wavelength').grid(row=y)
        wavelength = ttk.Entry(window, textvariable=self.wave).grid(row=y, column=1,ipadx=3, ipady=5, sticky='W')
        
        
        y+=1
        ttk.Label(window, text='distance').grid(row=y)
        wavelength = ttk.Entry(window, textvariable=self.pcdist).grid(row=y, column=1,ipadx=3, ipady=5, sticky='W')
        
        y+=1
        ttk.Label(window, text='S/N threshold').grid(row=y)
        wavelength = ttk.Entry(window, textvariable=self.sn_thresh).grid(row=y, column=1,ipadx=3, ipady=5, sticky='W')
        
        y+=1
        ttk.Label(window, text='Simulation? (T=1/F=0)').grid(row=y)
        wavelength = ttk.Entry(window, textvariable=self.sim).grid(row=y, column=1,ipadx=3, ipady=5, sticky='W')

        
        
        #temperature 
        y+=1
        ttk.Button(window, text='Get temp', command = lambda : self.load_temp_wcs()).grid(row=y, ipadx=5, ipady=15, sticky="ew")
        
        ttk.Entry(window, textvariable=self.temp, width=70).grid(row=y,column=1, ipadx=1, ipady=1,sticky='W')
        
        ttk.Button(window, text='Open fits', command = lambda : self.process_fits(self.temppath)).grid(row=y, column=2, ipadx=5, ipady=15, sticky="ew")
        
        # NH2
        y+=1
        ttk.Button(window, text='Get NH2', command = lambda : self.load_NH2()).grid(row=y, ipadx=5, ipady=15, sticky="ew")
        
        ttk.Entry(window, textvariable=self.nh2, width=70).grid(row=y,column=1, ipadx=1, ipady=1,sticky='W')
        
        ttk.Button(window, text='Open fits', command = lambda : self.process_fits(self.nh2path)).grid(row=y, column=2, ipadx=5, ipady=15, sticky="ew")
        
        ttk.Label(window, text='ImSeg Reduction Steps').grid(row=y, column =3, columnspan=2)

        
        # YSO
        y+=1
        ttk.Button(window, text='Get YSO', command = lambda : self.load_yso()).grid(row=y, ipadx=5, ipady=15, sticky="ew")
        
        ttk.Entry(window, textvariable=self.yso, width=70).grid(row=y,column=1, ipadx=1, ipady=1,sticky='W')
        
        ttk.Button(window, text='Open text file', command = lambda : self.read_yso()).grid(row=y, column=2, ipadx=5, ipady=15, sticky="ew")
        
        # RUN GET BLOBS
        self.getblobs = BooleanVar()
        self.getblobs.set(False)
        
        self.run_getblobs = ttk.Checkbutton(window, text="run Image\nSegmenation", var=self.getblobs).grid(row=y, column=3,sticky='W')
        
        # RUN NOISE BLOBS
        self.getnoiseblobs = BooleanVar()
        self.getnoiseblobs.set(False)
        self.run_getnoiseblobs = ttk.Checkbutton(window, text="run noise Image\nSegmentation", var=self.getnoiseblobs).grid(row=y, column=4,sticky='W')
        
        #YSO COLUMNS
        y+=1
        ttk.Label(window, text='YSO RA column #').grid(row=y)

        ttk.Entry(window, textvariable=self.racol).grid(row=y, column=1,ipadx=3, ipady=5, sticky='W')
        
        
        # RUN CORE SELECTION 
        self.coresel = BooleanVar()
        self.coresel.set(False)
        self.run_coreselection = ttk.Checkbutton(window, text="run core\nselection", var=self.coresel).grid(row=y, column=3, sticky='W', rowspan=2)
        
        
#        # RUN DO CORRECTIONS    
#        self.docorr = BooleanVar()
#        self.docorr.set(False)    
#        self.do_corrections = ttk.Checkbutton(window, text="do corrections", var=self.docorr).grid(row=y, column=4,sticky='W', rowspan=2)
        
        y+=1
        
        ttk.Label(window, text='YSO Dec column #').grid(row=y, column=0)
        
        ttk.Entry(window, textvariable=self.deccol).grid(row=y, column=1,ipadx=3, ipady=5, sticky='W')
        
        
#        y+=1
#        self.entry_var = StringVar()
#        ttk.Entry(window, textvariable=self.entry_var, width=20).grid(row=y,column=0, ipadx=1, ipady=1)
#        
#        ttk.Button(window, text='print noisepath / yso', command = lambda : self.process_filename(y)).grid(row=y, column=1,ipadx=5, ipady=15)
            
        # OUTBASE
        y+=1
        ttk.Button(window, text='Get outbase', command = lambda : self.get_outbase()).grid(row=y, ipadx=5, ipady=15, sticky="ew")
        
        ttk.Entry(window, textvariable=self.outbase, width=70).grid(row=y,column=1, ipadx=1, ipady=1,sticky='W')
        
        ttk.Button(window, text='Make outbase', command = lambda : self.make_outbase()).grid(row=y, column=2, ipadx=5, ipady=15, sticky="ew")
        
        y+=1
        ttk.Button(window, text='RUN', command = lambda : self.runImseg(), style='red.TButton').grid(row=y, column=3,ipadx=5, ipady=15)
        
        for row in range(y):
            window.grid_rowconfigure(row, weight=1)
        for col in range(3):
            window.grid_columnconfigure(col, weight=1)
            
        ttk.Separator(window, orient=VERTICAL).grid(column=2, row=0, rowspan=y, sticky='ns'+'E') 
        ttk.Separator(window, orient=VERTICAL).grid(column=3, row=0, rowspan=y, sticky='ns'+'W') 
        
    def load_signal(self):
        self.signalpath = askopenfilename(filetypes=(("FITS files", "*.fits"), ("All files", "*.*") ))
        self.signal.set(self.signalpath)
        
    def load_weight(self):
        self.weightpath = askopenfilename(filetypes=(("FITS files", "*.fits"), ("All files", "*.*") ))
        self.weight.set(self.weightpath)
        
    def load_sn(self):
        self.snpath = askopenfilename(filetypes=(("FITS files", "*.fits"), ("All files", "*.*") ))
        self.sn.set(self.snpath)
        
    def load_psf(self):
        self.psfpath = askopenfilename(filetypes=(("FITS files", "*.fits"), ("All files", "*.*") ))
        self.psf.set(self.psfpath)
        
    def load_noise_path(self):
        filename = askdirectory()
        self.noisepath.set(filename)
        
        
        
        
    def load_temp_wcs(self):
        temp = askopenfilename(filetypes=(("FITS files", "*.fits"), ("All files", "*.*") ))
        
        self.temppath = temp
        self.temp.set(self.temppath)
    
        self.wcs.set(self.temppath)
        
    def load_NH2(self):
        fil = askopenfilename(filetypes=(("FITS files", "*.fits"), ("All files", "*.*") ))
        self.nh2path = fil
        self.nh2.set(self.nh2path)
        
    def load_yso(self):
        fil = askopenfilename()
        self.ysopath = fil
        self.yso.set(self.ysopath)
        
    def read_yso(self):
        if self.ysopath != 'None': 
            print(np.loadtxt(self.ysopath))
            

    def get_outbase(self):
        filename = askdirectory()
        self.outbase.set(filename)
        
    def make_outbase(self):
        if not os.path.exists(self.outbase.get()):
            os.makedirs(self.outbase.get())
            print(f'{self.outbase.get()} folder is made')
        else:
            print('folder already exists')
        
        
    def process_fits(self, filename):
        os.system(f"ds9 {filename}")
        
    def process_filename(self,y):
        if self.entry_var.get() == 'noisepath':
            if self.noisepath.get():
                ttk.Label(window, text= self.noisepath.get()).grid(row=y, column=2)
            else:
                ttk.Label(window, text= "must provide path").grid(row=y, column=2)
        elif self.entry_var.get()  == 'yso':
            if self.yso.get():
                ttk.Label(window, text= self.yso.get()).grid(row=y, column=2)
            else:
                ttk.Label(window, text= "must provide file").grid(row=y, column=2)
        elif self.entry_var.get()  == 'outbase':
            if self.outbase.get():
                ttk.Label(window, text= self.outbase.get()).grid(row=y, column=2)
            else:
                ttk.Label(window, text= "must provide path").grid(row=y, column=2)
        elif self.entry_var.get()  == 'sim':
            if self.sim.get() is not None:
                if (self.sim.get() == 0) or (self.sim.get() == False):
                    ttk.Label(window, text= "False").grid(row=y, column=2)
                if (self.sim.get() == 1) or (self.sim.get() == True):
                    ttk.Label(window, text= "True").grid(row=y, column=2)
            else:
                ttk.Label(window, text= "must provide True or False").grid(row=y, column=2)
        else:
            ttk.Label(window, text= 'options are:\nnoisepath, yso, outbase, sim').grid(row=y, column=2)
            
    def runImseg(self):
        name = self.name.get()
        print('------------------- o ---------------------')
        print(name)
        signal = self.signal.get()
        sn = self.sn.get()
        weight = self.weight.get()
        psf = self.psf.get()
        noisepath = self.noisepath.get()
        wavelength = int(self.wave.get())
        pcdist = int(self.pcdist.get())
        sn_thresh = self.sn_thresh.get()
        sim = self.sim.get()
        temp = self.temp.get()
        wcs = self.wcs.get()
        coldens = self.nh2.get()
        ysoradec = self.yso.get()
        outbase = self.outbase.get()
        run_getblobs = self.getblobs.get()
        run_getnoiseblobs = self.getnoiseblobs.get()
        run_coreselection = self.coresel.get()
        ysora = self.racol.get()
        ysodec = self.deccol.get()
        instrument = self.instrument.get()

        
        os.system(f"python driver.py {name} {signal} {sn} {weight} {psf} {noisepath} {wavelength} {pcdist} {sn_thresh} {sim} {temp} {wcs} {coldens} {ysoradec} {outbase} {run_getblobs} {run_getnoiseblobs} {run_coreselection} {ysora} {ysodec} {instrument}")
    
        
        im = Image.open(f'{self.outbase.get()}/clumpfind_cores.png')
        
        width = 400
        height = 200

        im = im.resize((width,height), Image.ANTIALIAS)
        photo = ImageTk.PhotoImage(im)
        
        
        labels = ttk.Label(window, image = photo)
        labels.image = photo
        labels.grid(row=0, column=3,rowspan=4, columnspan=2)
        
        ttk.Button(window, text='Open blobprint fits', command = lambda : self.process_fits(f'{self.outbase.get()}/cores_blobprints.fits')).grid(row=4, column=3, ipadx=5, ipady=15, columnspan=2)
            
        if run_coreselection == True:
            im = Image.open(f'{self.outbase.get()}/histogram_cmf.png')
        
            width = 250
            height = 180

            im = im.resize((width,height), Image.ANTIALIAS)
            photo = ImageTk.PhotoImage(im)


            labels = ttk.Label(window, image = photo)
            labels.image = photo
            labels.grid(row=5, column=3,rowspan=7, columnspan=2)

    
    
        
window = Tk()
window.geometry("1300x700")
gui=Window(window)
window.mainloop() 
