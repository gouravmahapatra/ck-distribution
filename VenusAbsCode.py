#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 13:48:57 2018

Script to calculate F and P from a set of fourier files 
and a given set of geometry. 

@author: gouravmahapatr, TU Delft
"""
import os
import pymiedap.pymiedap as pmd
#import pymiedap.data_class as dt
import numpy as np
import matplotlib.pyplot as plt
#from scipy.signal import savgol_filter
from scipy import interpolate
import sys
# first load all the geometry information related to a particular orbit from the VEX data
path = '/Users/gouravmahapatr/Dropbox/PhD/spicav_data'
sys.path.append(path)
import read_spicav_ir as rd
sys.path.append('/Users/gouravmahapatr/Dropbox/PhD/Codes/PythonCodes')

#%%
"""
load a specific orbit on the basis of its orbit number 
"""

def find_nearest(array, value):
    ''' returns value,index'''
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx
    

def geom(index=-1):
    
    data = rd.processing_orbits(1478,8)
    # read the orbit geometries for the pass
    phase0 = data.geo.phase
    emission0 = data.geo.emission
    sza0 = data.geo.sza
    lat0 = data.geo.Lat
    lon0 = data.geo.Lon

    #phase = phase0[-2:-1]
    #emission = emission0[-2:-1]
    #sza = sza0[-2:-1]
    phase = phase0[np.array(index)]
    emission = emission0[np.array(index)]
    sza = sza0[np.array(index)]
    lat = lat0[np.array(index)]
    lon = lon0[np.array(index)]
    
    geo = [np.array([phase]),np.array([emission]),np.array([sza]),np.array([lat]),np.array([lon])]
    
    return geo


def flarr():
    # enter the starting lambda, d(lambda) and number of lambda
    lamb0 = 1.40050
    dlamb = 0.0005
    nlamb = 200
    # make sure the number of decimal places for file reading is correct 
    larr_str = [str('%1.5f'%(lamb0+dlamb*i)) for i in range(nlamb)]
    larr = np.array(larr_str,dtype=float)
    
    return larr_str,larr
    
    
def foukd2flux():    
    # get the geometry of the observation, 
    # default geometry is closest to equator
    phase,emi,sza,lat,lon = geom(index=0)
    
    # point the function to the fourier folder location
    set_path0 = '/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/fullKDIS/DAPvenusCombLay/FourierFiles'
    os.chdir(set_path0)
    
    larr,nlarr = flarr()
    
    I,Q,U,V = [],[],[],[]
    for j in range(len(larr)):
        print('Wavelength: ',larr[j])
        # enter the folder with the wavelength's name
        set_path1 = set_path0+'/'+larr[j]
        os.chdir(set_path1)
        
        Il,Ql,Ul,Vl = [],[],[],[]
        
        for i in range(10):
            fname = 'fou'+str(i+1)+'.out'
            #print('Now reading... ',fname)
            
            # use pymiedap's read_dap_output to get the flux vector
            I0,Q0,U0,V0 = pmd.read_dap_output(phase,sza,emi,fname)
            
            # append these values to the stokes vectors lists
            # multiply the stokes vectors with the Gaussian weights. 
            Il.append(I0[0]*kdgw[1,i])
            Ql.append(Q0[0]*kdgw[1,i])
            Ul.append(U0[0]*kdgw[1,i])
            Vl.append(V0[0]*kdgw[1,i])
        
        I.append(sum(Il))
        Q.append(sum(Ql))
        U.append(sum(Ul))
        V.append(sum(Vl))
        
    # convert to array
    I = np.array(I)
    Q = np.array(Q)
    U = np.array(U)
    V = np.array(V)
    
    return I,Q,U,V
    
    
def getSPIRflux(lmin=None, lmax=None,idx=0,det=0):
    """
    This function gets the NORMALIZED flux vectors from the SPICAV-IR routine
    within a requested given range of wavelengths.
    You can input an orbital epoch also, default is set at index=0.
    
    lmin and lmax are in nanometers!!
    """
    data = rd.processing_orbits(1478,8)
    
    if det == 0:
        r0 = data.r0[:,idx]
        w0 = data.w0[:,idx]
    else:
        r0 = data.r1[:,idx]
        w0 = data.w1[:,idx]
        
    # drop the "nan's" and normalize the radiance
    w0 = w0[~np.isnan(r0)]
    r0 = r0[~np.isnan(r0)]  # remove nan's
        
    if (lmin==None) & (lmax == None):
        # get the SPICAV-IR radiance and wavelength data
        r = r0
        w = w0
    else:
        r = r0[(w0>lmin)&(w0<lmax)]
        w = w0[(w0>lmin)&(w0<lmax)]
        
    r = r/np.max(r)
        
    return r,w
    
def getSPIRpol(lmin=None, lmax=None,idx=0,det=0,orbn=1478,orba=8):  
    """
    This function returns polarized flux data for a given orbit geometry =idx.
    
    It does not do any processing as done by Loic Rossi to arrive at his selected 
    bands.
    
    pol = (det1 - det0)/(det0 + det1)
    
    Returns: pol,wavelength
    
    Authr: G. Mahapatra
    
    """
    
    data = rd.processing_orbits(orbn,orba)
    
    r0 = data.r0[:,idx]
    w0 = data.w0[:,idx]
    
    r1 = data.r1[:,idx]
    w1 = data.w1[:,idx]
    
    
    if (lmin==None) & (lmax == None):
    # get the SPICAV-IR radiance and wavelength data
        r0 = r0
        w0 = w0
        r1 = r1
        w1 = w1
    else:
        lidx = (w0>lmin)&(w0<lmax)
        r0 = r0[lidx]
        w0 = w0[lidx]
        r1 = r1[lidx]
        w1 = w1[lidx]
        
    # calculate the polarization at all the available geometries
    # polarization is calculated as (det1 - det0)/(det0 + det1)
    pol = np.nan_to_num((r1 - r0)/(r0 + r1))
    
    return pol,w1

 #%%      
def spicav_irf(perch):
    ''' This function provides the SPICAV IRF as a function of wavenumber
    interval. '''
    # load the SPICAV IR PSF file
    path = '/Users/gouravmahapatr/Dropbox/PhD/spicav_data/psf_lw_all_O2_1270_desc.txt'
    irf_data = np.loadtxt(path,skiprows=1)
    dnu = irf_data[:,0]             # in cm-1
    df = irf_data[:,1]              # in kHz
    ch0 = irf_data[:,2]             # this is the normalized instrument response
    
    ch0_2 = np.array(ch0)
    ch0_2[(dnu>-4.5)&(dnu<6)] = 0
    
    #perch = -20.0   # this is the percentage change normalized to 1.
    ch0_3 = ch0 + perch*ch0_2/1e2
    #plt.plot(dnu,ch0_3,label=perch)
    
    # make a new range of dnu
    dnu_new = np.arange(-100,100,abs(dnu[0]-dnu[1]))
    
    # make the interpolator to get the instruments response for a given wavenumber 
    ch0_intp = interpolate.Akima1DInterpolator(dnu,ch0_3)
    
    ch0_4 = ch0_intp(dnu_new)
    
    # convert the nans into zeros
    ch0_4[np.isnan(ch0_4)] = 0
    
    # make a new interpolation box
    ch0_intp = interpolate.Akima1DInterpolator(dnu_new,ch0_4)
    
    # these are the constants to determine the wavenumber (Korablev et al.,)
    f = 1.e3
    a0 = -3.3865473e-8
    b0 = 7.2595705e-2
    c0 = -2.0449838e+0
    a1 = -3.5371703e-8
    b1 = 7.2919764e-2
    c1 = -1.9140569e+1     
    # given a frequency of operation determine wavenumber
    nu_0 = a0*f**2 + b0*f + c0   # in cm-1
    nu_1 = a1*f**2 + b1*f + c1   # in cm-1
    # given a wavenumber determine the frequency 
    
    # model the response function as (sin(x)/x)^2.
#    x = dnu
#    f = 14 # in cm-1
#    xd = (2*np.pi*x/f)+np.pi*1.2
#    y = (np.sin(xd)/(xd*1.0))**2
#    plt.plot(x,y)
    return ch0_intp,ch0_3,dnu_new,dnu,ch0
#%%

def wvl2wvn(wavelength):
    ''' This function converts wavelength (in microns!!!) to 
    wavenumber (in cm-1).                              '''
    wavenumber = 1.e4/wavelength
    return wavenumber

def wvn2wvl(wavenumber):
    ''' This function converts wavenumber (in cm-1) 
    to wavelength (in microns!!!).              '''
    wavelength = 1.e4/wavenumber
    return wavelength

#%%
def slitfunction(wavel0,sigma,nw,truncw):
    ''' This function determines the wavenumber eqivalent of
    wavelength values and constructs a slitfunction with a given 
    width (sigma).
    
    It outputs min and max wavenumbers (vmin & vmax), and 
    truncation limits (tmin & tmax) depending on the truncation 
    wavelength (truncw).
    
    specv: the edges of the domain in wavenumber
    
    '''
    wavelmin = wavel0 - 0.5*sigma
    wavelmax = wavel0 + 0.5*sigma
    
    if wavelmin < 0.0: 
        raise ValueError('Error: wavelmin out of bounds!')
    
    # number of wavenumber values (nnv)
    nnv = int((wavelmax-wavelmin)*nw)
    
    # convert the wavelength into wavenumber 
    vmin = wvl2wvn(wavelmax) # maximum wavel goes to min wavenumber!
    vmax = wvl2wvn(wavelmin)  
    
    # truncation in wavenumbers (max. value is chosen)
    truncv = -1.e4*(1.0/(wavelmax+truncw)-1.0/wavelmax)
    if truncv < 0.0: truncv=1.e4/wavelmax
    
    # these are truncation variables that are used in reading the line param file!
    tmin = vmin - truncv
    tmax = vmax + truncv
    
    # fill the specv array with equally spaced wavenumbers
    specv = []
    i,specvi  = 0,0
    while specvi < vmax:
        specvi = vmin+(vmax-vmin)*i/nnv
        specv.append(specvi)
        i+=1
    nnv = i
    specv = np.array(specv,dtype=float)
    specw = wvn2wvl(specv)
    
    # fill the speci array which is later used in quicksort array
    # speci stores the integrated slitfunction.
    # this really looks like a normalized array. Need to look further. 
    tot = 0
    speci = np.zeros(nnv)
    for i in range(nnv-1,0,-1):
        w1 = specw[i-1]
        w2 = specw[i]
        plus = (w1-w2)
        tot = tot+plus
        speci[i]=tot
    speci = speci/tot
        
    return nnv,truncv,specv,speci,vmin,vmax,tmin,tmax

#%%

def convolve_spectra(spectra,wav,wav_conv,perch):    
    ''' Find the integrated spectral point at a given wavelength value.
        This function convolves the line-by-line spectra with SPICAV-IR.
        
        As of now it's coded to take an interval specified by 'dnu' of 
        the instrument. The instruments response changes rapidly with wavelength
        so this needs to be fixed. 
        
        Inputs: spectra = as defined for the line-by-line high res
                wav = wavelength arra of the high res 
                wav_conv = wavelength points to be convolved around.
                perch = percentage change of the instrument response             
    ''' 
    # get the equivalent wavenumbers
    nu = 1e4/wav
    nu_conv = 1e4/wav_conv
    
    # load the instrument response function
    ch0_intp,_,_,dnu,_ = spicav_irf(perch)
    #irf = gaussian_IRF(specv,sig=sig,mu=waven0)
    
    # fix the number of elements to be taken around the central wavenumber
    nelm = len(nu[nu <= (nu[-1]+100)])
    
    I = spectra[:,:,1]
    Q = spectra[:,:,2]
    
    # separate the initial I and Q to form the angular fluxes.
    I_0 = 0.5*(I+Q)
    I_90 = 0.5*(I-Q)
    
    #Ich = []
    Ich_0 = []
    Ich_90 = []
    c=0
    
    length_test = []
    
    for i in range(len(nu)):
        # determine the central wavenumber (nu0)
        nu0 = nu[i]
        
        # find the nearest wavenumber nu and index to nu0
        _,idx_nu0 = find_nearest(nu,nu0)
        
        # find the min and max in wavenumber on the basis of instruments d(nu).
        nu_min = nu0+dnu.min()
        nu_max = nu0+dnu.max()
        # corresponding wavelengths
        wmin = 1e4/nu_max
        wmax = 1e4/nu_min
        
        # make an equally spaced array with resolution specified in terms of nu
        #nu_slit = nu[(nu>=nu_min)&(nu<nu_max)]
        left_nelm = i - nelm
        right_nelm = i + nelm
        if i < nelm:
            left_nelm = 0
            right_nelm = i
        if i == 0:
            left_nelm = 0
            right_nelm = 1
        if right_nelm >= (len(nu)-1):
            right_nelm = (len(nu)-1)
            left_nelm = i - ((len(nu)-1) - i)
        if left_nelm <= 0:
            left_nelm = 0
            right_nelm = i+i
        if i == (len(nu)-1):
            right_nelm = (len(nu)-1)
            left_nelm = (len(nu)-2)
        
        # to handle the edge-cases
        if i == 0:
            right_nelm = 2
        elif i == (len(nu)-1):
            left_nelm = -3
            right_nelm = -1
            
        nu_slit = nu[left_nelm:right_nelm]

        # find the values in the spectra that lie within the specv range
        I_0_slit = I_0[left_nelm:right_nelm,:]
        #plt.scatter((1e4/nu0),len(I_0_slit),s=2)
        I_90_slit = I_90[left_nelm:right_nelm,:]
                            
        #dnu = abs(nu_slit[0]-nu_slit[1]) # in wavenumbers
        dnu_slit = nu_slit - nu0
        irf_slit = ch0_intp(dnu_slit)
        irf_slit[np.isnan(irf_slit)] = 0
        
        #I_slit = spectra[(nu>=nu_min)&(nu<nu_max),0,1]
        
        # find the flux values inside the bins of the slit and sum them up
        #I_slit_afterirf = I_slit*irf_slit
        
        # make an irf box with shape equal to the Fluxes
        irf_slit_box = np.zeros(I_0_slit.shape)
        #print(I_0_slit.shape)
        for j in range(len(I_0_slit[0,:])):
            irf_slit_box[:,j] = irf_slit
        
        I_0_slit_afterirf = I_0_slit*irf_slit_box
        I_90_slit_afterirf = I_90_slit*irf_slit_box
        
        #length_test.append(len(I_0_slit_afterirf))
        
        Ich_0.append(np.sum(I_0_slit_afterirf,axis=0))
        Ich_90.append(np.sum(I_90_slit_afterirf,axis=0))
        
        if c == 3000:
            print('Wavelength process number is ',i)
            c=0
        
        c+=1
        
    # recombine to form I and Q
    Iconv = np.array(Ich_0)+np.array(Ich_90)
    Qconv = np.array(Ich_0)-np.array(Ich_90)
    
    # pack it into a spectra array
#    spectra = np.zeros((Iconv.shape[0],Iconv.shape[1],3))
#    spectra[:,:,1] = Iconv
#    spectra[:,:,2] = Qconv
    
    return Iconv,Qconv

#%%
def convolve_spectra(spectra,wav,wav_conv,perch):    
    ''' Find the integrated spectral point at a given wavelength value.
        This function convolves the line-by-line spectra with SPICAV-IR.
        
        As of now it's coded to take an interval specified by 'dnu' of 
        the instrument. The instruments response changes rapidly with wavelength
        so this needs to be fixed. 
        
        Inputs: spectra = as defined for the line-by-line high res
                wav = wavelength arra of the high res 
                wav_conv = wavelength points to be convolved around.
                perch = percentage change of the instrument response             
    ''' 
    # get the equivalent wavenumbers
    nu = 1e4/wav
    nu_conv = 1e4/wav_conv
    
    # load the instrument response function
    ch0_intp,_,_,dnu,_ = spicav_irf(perch)
    #irf = gaussian_IRF(specv,sig=sig,mu=waven0)
    
    I = spectra[:,:,1]
    Q = spectra[:,:,2]
    spectra_conv = np.zeros(spectra.shape)
    
    # initialize the convolution arrays
    Iconv = np.zeros((len(nu_conv),I.shape[1]))
    Qconv = np.zeros((len(nu_conv),I.shape[1]))
    
    for i in range(len(nu_conv)):
        I_temp = []
        Q_temp = []
        
        nu0 = nu_conv[i]
        print(1e4/nu0)
        
        # find the min and max in wavenumber on the basis of instruments d(nu).
        nu_min = nu0+dnu.min()
        nu_max = nu0+dnu.max()
        # corresponding wavelengths
        wmin = 1e4/nu_max
        wmax = 1e4/nu_min
        
        # find all the values in the line-by-line lying within the nu limit
        idx_nu_lim = np.array(np.where((nu>nu_min)&(nu<nu_max)))
        nu_lim = nu[np.squeeze(idx_nu_lim)]
        
        # get the dnu values from nu
        dnu_from_nu = np.array((nu_lim - nu0))
        
        # get the irf frm the dnu values
        irf = ch0_intp((dnu_from_nu))
        
        for j in range(I.shape[1]):  # this loop is across the geometries
            I_temp.append(np.average(irf*np.squeeze(I[idx_nu_lim,j])))
            Q_temp.append(np.average(irf*np.squeeze(Q[idx_nu_lim,j])))
            
        I_temp = np.array(I_temp)
        Q_temp = np.array(Q_temp)
        
        plt.scatter(1e4/nu0,I_temp[30])
        
        Iconv[i,:] = I_temp
        Qconv[i,:] = Q_temp
        
    # pack it into a spectra array
    spectra_conv = np.zeros((Iconv.shape[0],Iconv.shape[1],3))
    spectra_conv[:,:,1] = Iconv
    spectra_conv[:,:,2] = Qconv
    
    return spectra_conv  
        
#%%

def kd_spectrum(filename=None,wav_file=None,atm_file=None):
    """
    load the spectrum data used for making k-distribution.
    This spectrum has a different arrangement as compared
    to a continuous spectra. 
    
    Shape is determined as: [nwav,nnv,nlev] where,
    nwav = no. of central wavelengths
    """
    if filename is None and wav_file is None and atm_file is None:
        filename='/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/fullKDIS/KDIS_OK/spectrum.out'
        wav_file='/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/fullKDIS/KDIS_OK/wav.dat'
        atm_file='/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/fullKDIS/KDIS_OK/VenusIg.dat'
        kd_file='/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/fullKDIS/KDIS_OK/kd.out'
        
    #load the levels data from 'venusIg.dat' file
    atm_data = [line.strip().split() for line in open(atm_file)]  # set the atmosphere model
    nlev = np.array([int(atm_data[i][0]) for i in range(6,len(atm_data))])
    plev = [np.float(atm_data[i][1]) for i in range(6,len(atm_data))] # in bars
    tlev = [np.float(atm_data[i][2]) for i in range(6,len(atm_data))]
        
    wavdat = np.loadtxt(wav_file)  # wavelengths used in microns
    
    spec_data = [line.strip().split() for line in open(filename)]
     
    kdata = np.loadtxt(kd_file)
    
    nnv = int(spec_data[8+len(nlev)][3])   # number of intervals per central wavelength
    
    ngp = 20   # number of Gauss points in the x-interval!!
    
    if np.ndim(wavdat) == 1:
        wav = wavdat[1]
        spec = np.zeros((1,nnv,len(nlev)))
        abs_spec = np.zeros((1,nnv,len(nlev)))
        dim_wav = 1
    else:
        wav = wavdat[:,1]
        spec = np.zeros((len(wav),nnv,len(nlev)))
        abs_spec = np.zeros((len(wav),nnv,len(nlev)))
        dim_wav = len(wav)
        kdis = np.zeros((len(nlev),len(wav),ngp))
 
    k = 10+len(nlev)
    
    kcount = 0
    for l in range(int(len(kdata)/dim_wav)):
        for i in range(dim_wav):          
            for j in range(nnv):
                if j == (nnv-1):
                    spec[i,j,l] = np.float128(spec_data[k][1])
                    abs_spec[i,j,l] = np.float128(spec_data[k][2])
                    k = k+5
                else:
                    k = k+1
                    spec[i,j,l] = np.float128(spec_data[k][1])
                    abs_spec[i,j,l] = np.float128(spec_data[k][2])
            kdis[l,i,:] = kdata[kcount,4:]
            kcount+=1
            
    return spec,abs_spec,kdis


def rotate_stokes(Q,U,beta):
    '''This function takes in the Stokes vector elements Q,U
    and performs a rotation by an angle beta which is the 
    angle between the planetary scattering plane and local scattering
    plane.
    Input:
        Q,U: Stokes elements (can be arrays)
        beta [degrees!]
    '''
    beta_r = np.deg2rad(beta)           # convert beta into radians 
    Qnew = np.cos(2*beta_r)*Q+np.sin(2*beta_r)*U
    Unew = -np.sin(2*beta_r)*Q+np.cos(2*beta_r)*U
    
    return Qnew,Unew
#%%
def gaussian_IRFlbl(x,sig=1.0,mu=0):
    ''' 
    In VenusAbsCode.py:
        
    Given a list of x-values, calculate Gaussian dist. 
    "sig" can be used to specify the std.dev. of the function in nm.
    '''
    #sigv = abs(wvl2wvn(1.000)-wvl2wvn(1.000-sig*1e-3))
    sigv2 = sig**2
    mu = 0.0
    
    p_x = (1/np.sqrt(2*np.pi*sigv2))*np.exp(-((x-mu)**2/(2*sigv2)))
    return p_x

# %%
def convolve_spectra_gaussian(spectra,wav,wav_conv,var=0.5,irf='box',sigma=1.e-3):    
    ''' Find the integrated spectral point at a given wavelength value.
        This function convolves the line-by-line spectra with a model Gaussian.
        
        The Gaussian's response width is fixed for now. 
        
        If irf = 'box' (default for now) then the Gaussian response is surpassed.
        
        Inputs:
            spectra: A specific type of constructed array containing line-by-line
                            Stokes parameters at given geometries and wavelengths.
            wav:    list of wavelengths used in line-by-line!!!
            var:    the variance of the Gaussian IRF
            
    '''    
    # set the slit function width (in microns)
    # sigma = 1*1.e-3
      
    # set the gaussian width
    sig = var*(20/sigma)

    # determine the approximate number of elements in the high-res array
    # corresponding to this sigma value
    nelm = len(wav[wav<=(wav[0]+sigma)])
    
    # set the resolution per micron
    nw = 500000
    
    # set the truncation window (in microns)
    truncw = 0.2 
    
    # get the equivalent wavenumbers
    nu = 1e4/wav
        
    I = spectra[:,:,1]
    Q = spectra[:,:,2]
    U = spectra[:,:,3]
    
    
    Iconv = []
    Qconv = []
    Uconv = []
    
    c=0
    
    for i in range(len(wav_conv)):
        # get the central wavelength
        wavel0 = wav_conv[i]
        # determine the central wavenumber (nu0)
        waven0 = nu[i]
        
        # get the index and nearest wavelength in the high-res file
        w0,idx_w0 = find_nearest(wav,wavel0-0.5*sigma)
        
        # get the slit characteristics
        nnv,truncv,specv,speci,vmin,vmax,tmin,tmax = slitfunction(wavel0,sigma,nw,truncw)
        specw = wvn2wvl(specv)
        wmin = min(specw)
        wmax = max(specw)
        
        #irf = np.ones(irf.shape)
        # make an interpolation box
        #irf_intp = interpolate.interp1d(specv,irf,kind='cubic')
        #plt.plot(irf)           
        # make an equally spaced array with resolution specified in terms of nu
        #wav_slit = wav[(wav>=wmin)&(wav<=wmax)]      
#        irf_slit = irf_intp(nu_slit)
            
        # find the values in the spectra that lie within the specv range
        idxf = (idx_w0+nelm)         # this gets the ending index number
        if idxf >= (len(I)-1):
            idxf = len(I)-1
        
        I_slit = I[idx_w0:idxf,:]
        Q_slit = Q[idx_w0:idxf,:]
        U_slit = U[idx_w0:idxf,:]
        
        # load the instrument response function
        # this simulates a box function
        if irf == 'box':
            irf_slit = np.ones(I_slit.shape[0])
            print('box IRF activated!')
        else:
            x = np.linspace(-10,10,len(I_slit))
            irf = gaussian_IRFlbl(x,sig=sig,mu=0.0) 
            # normalize the irf
            irf_slit = irf/max(irf)
        
        # find the flux values inside the bins of the slit and sum them up
        I_slit_afterirf = np.transpose(I_slit)*irf_slit
        #plt.plot(wav_slit,I_slit_afterirf,ls='dashed')
        #plt.scatter(wavel0,len(I_slit_afterirf),s=2.5,c='k')
        Q_slit_afterirf = np.transpose(Q_slit)*irf_slit
        U_slit_afterirf = np.transpose(U_slit)*irf_slit

        # sum/average them all to get one pixel value for each wavelength value
        ndata = len(I_slit_afterirf)
        Iconv.append(np.average(I_slit_afterirf,axis=1))
        Qconv.append(np.average(Q_slit_afterirf,axis=1))
        Uconv.append(np.average(U_slit_afterirf,axis=1))
        
        if c == 3000:
            print('Wavelength process number is ',i)
            c=0
        
        c+=1
        
    # convert into array
    Iconv = np.array(Iconv)
    Qconv = np.array(Qconv)
    Uconv = np.array(Uconv)
    
    # create the convolved spectra
    spectra_convolved = np.zeros([len(wav_conv),spectra.shape[1],spectra.shape[2]])
    spectra_convolved[:,:,1] = Iconv#/np.max(Iconv,axis=0)
    spectra_convolved[:,:,2] = Qconv#/np.max(Qconv,axis=0)
    spectra_convolved[:,:,3] = Uconv#/np.max(Uconv,axis=0)
    spectra_convolved[:,:,4] = np.sqrt(spectra_convolved[:,:,2]**2+spectra_convolved[:,:,3]**2)/spectra_convolved[:,:,1]
    
    print('convolve_spectra_gaussian: The line-by-line spectra is normalized!!!')
#    return Iconv,Qconv,Uconv    
    return spectra_convolved

#%%
def convolve_spectra_box(spectra,wav,wav_conv,sigma=2.e-3):
    '''
    This function specifically convolves the high-resolution line-by-line 
    spectra with a 'box' convolution function.
    sigma:     window around the central wavelength (in microns!!!)
    wav:       wavelength list of the line-by-line spectra. Dimension should match 
                the spectra wavelength dimension.
    wav_conv:  wavelength list containing central wavelengths around which
               the convolution should be made.
    '''
    
    # get the Stokes vectors from the spectra array
    I = spectra[:,:,1]
    Q = spectra[:,:,2]
    U = spectra[:,:,3]
    
    # arrays where the convolved values shall be stored
    Iconv = []
    Qconv = []
    Uconv = []
    
    # file counter
    c = 0
    
    for i in range(len(wav_conv)):
        # get the central wavelength
        wavel0 = wav_conv[i]
        
        # get the index of the values in line-by-line within the limits of the 
        # central wavelength
        idx_slit = np.squeeze(np.where((wav>=(wavel0-sigma*0.5))&(wav<=(wavel0+sigma*0.5))))
    
        # check if the idx_slit got any values. If not, then add nan's to the array
        if len(idx_slit) == 0:
            print('Wavelength value ',wavel0,' does not have any line-by-line spectra!')
        
        # Stokes values in the slit
        I_slit = I[idx_slit,:]
        Q_slit = Q[idx_slit,:]
        U_slit = U[idx_slit,:]
        
        # average them across the wavelength interval
        Iconv.append(np.mean(I_slit,axis=0))
        Qconv.append(np.mean(Q_slit,axis=0))
        Uconv.append(np.mean(U_slit,axis=0))
        
        # print progress...
        if c == 3000:
            print('Wavelength process number is ',i)
            c=0   
        c+=1
    
    # convert into array
    Iconv = np.array(Iconv)
    Qconv = np.array(Qconv)
    Uconv = np.array(Uconv)
    
    # create the convolved spectra
    spectra_convolved = np.zeros([len(wav_conv),spectra.shape[1],spectra.shape[2]])
    spectra_convolved[:,:,1] = Iconv#/np.max(Iconv,axis=0)
    spectra_convolved[:,:,2] = Qconv#/np.max(Qconv,axis=0)
    spectra_convolved[:,:,3] = Uconv#/np.max(Uconv,axis=0)
    spectra_convolved[:,:,4] = (-spectra_convolved[:,:,2])/spectra_convolved[:,:,1]
    
    return spectra_convolved  


#%% 
def convolve_solarflux(wavin):
    '''
    This function returns the solar spectra at the exact points as supplied by a wavelength. 
    Wavelength should be in microns!
    '''
    # get the path to the file
    solarflux = np.loadtxt('/Users/gouravmahapatr/Dropbox/PhD/Codes/PythonCodes/solarflux_venus.dat')
    # make an interpolation object as a function of solar flux and wavelength
    wavsun = solarflux[:,0]*1.e-3       # converted into microns!!!
    flux = solarflux[:,1]
    intp_solar = interpolate.interp1d(wavsun,flux)   
    return intp_solar(wavin)
#%%
# %%

