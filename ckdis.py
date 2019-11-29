#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 12:52:10 2018

KDIS python code, VERIFIED!!

This code can be used to generate the absorption cross-sections of CO$_2$ and 
then perform a correlated k-distribution of the absorption cross-sections 
to simulate the signal in a given wavelength range.

Original code (in FORTRAN) written by Dr. D.M. Stam, TU Delft.

@author: gouravmahapatr, TU Delft
"""

import numpy as np
import matplotlib.pyplot as plt 
from VenusAbsCode import *
import os
import sys
from scipy import interpolate
import linestyles as ls
import time
#%% This section decides on which line files to use
'''
We have two options, 
    1) is locally downloaded '.par' file. Option: 'local'
    2) retrieve the latest HITRAN line files from the server. Option: 'hapi'

For option 2), you need to have the 'hapi.py' program in your working path.
'hapi.py' link: https://hitran.org/static/hapi/hapi.py(last checked, 29.11.2019)
    
'''
if 'tableList' in globals():
    if 'CO2' not in tableList():
        from hapi import *
elif 'tableList' not in globals():
    from hapi import *

kdis_path = '/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/fullKDIS/KDIS_OK/'
os.chdir(kdis_path)
sys.path.append(kdis_path)

import salpha
import spectrum

'''
read the whole line param file at once 
'''
# this variable specifies which source to read from
# 'hapi' = HITRAN python reader
# 'local' = locally written script
line_source = 'local' 
#line_source = 'hapi'
print('The',line_source,' line source is being used to read the line parameters!!!')

atm_profile = np.loadtxt(kdis_path+'VenusIg.dat',skiprows=6)
line_fname = 'CO2_hit04.par' #'5c181ff0.par'
#line_fname = '5d137619.par'
#line_fname = '5d331caa.par'
global line_data
if 'line_data' not in globals() and line_source == 'local':
    print('Reading absorption parameters from ',line_fname,'.')
    file_name = line_fname
    line_data = [line for line in open(kdis_path+file_name)]

# this is where the file from HITRAN server is downloaded.
if os.path.isfile('CO2.data')==True and line_source == 'hapi':
    if 'tableList' in globals():
        if 'CO2' not in tableList():
            fetch_by_ids('CO2',[7,8,9,10,11,12,13,14,121,15,120,122],200,50000)
elif os.path.isfile('CO2.data')==False and line_source == 'hapi':
    fetch_by_ids('CO2',[7,8,9,10,11,12,13,14,121,15,120,122],200,50000)


#%% This section has the 'optional' .par file format for reference,
#   and also the line ID and some declarations. 

istout = 4
ngas2imol = 2  # this is the gas ID. Needs to be incorporated in future. CO2 = 2.

# this is for the fortran type format. taken from max_incl file.
nvMAX=661000
nlevMAX=40
ngaussMAX=50
ndimMAX=200
nwavsMAX=5000

# converts between TRANSPORT_FORMAT and OBJECT_FORMAT
HITRAN_FORMAT_160 = {
   'M'          : {'pos' :   1,   'len' :  2,   'format' : '%2d' },
   'I'          : {'pos' :   3,   'len' :  1,   'format' : '%1d' },
   'nu'         : {'pos' :   4,   'len' : 12,   'format' : '%12f'},
   'S'          : {'pos' :  16,   'len' : 10,   'format' : '%10f'},
   'R'          : {'pos' :  26,   'len' :  0,   'format' : '%0f' },
   'A'          : {'pos' :  26,   'len' : 10,   'format' : '%10f'},
   'gamma_air'  : {'pos' :  36,   'len' :  5,   'format' : '%5f' },
   'gamma_self' : {'pos' :  41,   'len' :  5,   'format' : '%5f' },
   'E_'         : {'pos' :  46,   'len' : 10,   'format' : '%10f'},
   'n_air'      : {'pos' :  56,   'len' :  4,   'format' : '%4f' },
   'delta_air'  : {'pos' :  60,   'len' :  8,   'format' : '%8f' },
   'V'          : {'pos' :  68,   'len' : 15,   'format' : '%15s'},
   'V_'         : {'pos' :  83,   'len' : 15,   'format' : '%15s'},
   'Q'          : {'pos' :  98,   'len' : 15,   'format' : '%15s'},
   'Q_'         : {'pos' : 113,   'len' : 15,   'format' : '%15s'},
   'Ierr'       : {'pos' : 128,   'len' :  6,   'format' : '%6s' },
   'Iref'       : {'pos' : 134,   'len' : 12,   'format' : '%12s'},
   'flag'       : {'pos' : 146,   'len' :  1,   'format' : '%1s' },
   'g'          : {'pos' : 147,   'len' :  7,   'format' : '%7f' },
   'g_'         : {'pos' : 154,   'len' :  7,   'format' : '%7f' }
}

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
def rdabsfile(vdmin=None,vdmax=None,trunc=0):
    ''' This function reads the HITRAN 2004 line transitions file...
    It returns 9 parameters back as a function of wavenumber vd (cm-1).
    i2 = molecule ID
    id = isotopologue number 
    vd = frequency of transition in wavenumbers (cm-1)
    rsd = intensity in cm-1/(molecule cm-2) at 296 K
    wtm = Einstein A-coefficient in s-1.
    rad = air-broadenend halfwidth in cm-1/atm at 296 K
    sbh = self-broadened Lorentzian HWHM in cm-1/atm
    red = lower state energy in wavenumbers (cm-1)
    rtd = coefficient of temperature dependence
    
    ntot = number of transition lines within [tmin,tmax].
    
    It is an adaptation from the original code written by D.M. Stam (TU Delft).
    '''
    # initialization
    i2,ID,vd,rsd,wtm,rad,sbh,red,rtd = [],[],[],[],[],[],[],[],[]
    for i in range(1,len(line_data)):
        i2.append(int(line_data[i][0:2]))
        # for isotopologue ID
        if line_data[i][2:3] == 'A':    
            ID.append(11)  # 'A' = isotopologue ID = 11
        elif line_data[i][2:3] == 'B':
            ID.append(12)  # 'B' = isotopologue ID = 12
        else:
            ID.append(int(line_data[i][2:3]))
        #ID.append(int(line_data[i][2:3]))
        vd.append(float(line_data[i][3:15]))
        rsd.append(float(line_data[i][15:25]))
        wtm.append(float(line_data[i][25:35]))
        rad.append(float(line_data[i][35:40]))
        sbh.append(float(line_data[i][40:45]))
        red.append(float(line_data[i][45:55]))
        rtd.append(float(line_data[i][55:59]))
    
    # convert into arrays    
    i2 = np.array(i2)
    ID = np.array(ID)
    vd = np.array(vd)
    rsd = np.array(rsd,dtype=float)
    wtm = np.array(wtm)
    rad = np.array(rad)
    sbh = np.array(sbh)
    red = np.array(red)
    rtd = np.array(rtd)
    
    if (vdmin==None) & (vdmax==None):
        print('Warning: reading the whole HITRAN file as no limits provided!')
    else:
        # determine the limits for the wavenumber cutoff
        tmin = vdmin - trunc
        tmax = vdmax + trunc
        if (tmin < min(vd)) or (tmax > max(vd)):
            raise ValueError('truncation exceeds available line parameters extent.')
    
        i2 = i2[(vd>tmin)&(vd<tmax)]
        ID = ID[(vd>tmin)&(vd<tmax)]
        vd = vd[(vd>tmin)&(vd<tmax)]
        rsd = rsd[(vd>tmin)&(vd<tmax)]
        wtm = wtm[(vd>tmin)&(vd<tmax)]
        rad = rad[(vd>tmin)&(vd<tmax)]
        sbh = sbh[(vd>tmin)&(vd<tmax)]
        red = red[(vd>tmin)&(vd<tmax)]
        rtd = rtd[(vd>tmin)&(vd<tmax)]
   
        
    ntot = len(i2) # this variable passes the number of lines in the requested range
        
    return i2,ID,vd,rsd,wtm,rad,sbh,red,rtd,ntot
#%%
def rdHAPIparam():
    '''
    This function downloads and reads the latest HITRAN line files using the 
    HAPI code.
    
    Ref:
    '''
    # import only the essential functions for now
    #
    filename = 'CO2'
    # this command fetches all the listed CO2 isotopologues listed 
    # on the HITRAN website and stores them in the current directory. 
    # The last two arguments specify the desired wavenumber range. 
    # for help type getHelp(fetch_by_ids)

    # this command describes the downloaded table.
    describeTable('CO2')
    # retrieve all the 9 parameters that are also in HITRAN 2004.
    molec_id = getColumn('CO2','molec_id')
    local_iso_id = getColumn('CO2','local_iso_id')
    nu = getColumn('CO2','nu')
    sw = getColumn('CO2','sw')
    a = getColumn('CO2','a')
    gamma_air = getColumn('CO2','gamma_air')
    gamma_self = getColumn('CO2','gamma_self')
    elower = getColumn('CO2','elower')
    n_air = getColumn('CO2','n_air')
    
    # now equate these HAPI variables to the old Stam style variables.
    i2 = molec_id
    ID = local_iso_id
    vd = nu
    rsd = sw
    wtm = a
    rad = gamma_air
    sbh = gamma_self
    red = elower
    rtd = n_air
    
    # total number of lines found
    ntot = len(rsd) 
    
    return i2,ID,vd,rsd,wtm,rad,sbh,red,rtd,ntot
    
#%%
def abs_lim_check(vdmin,vdmax,truncv,vd):
    '''
    This function checks if the central wavelength lies within the limit 
            [wavel0 - trunc],[wavel0+trunc].
    '''
    tmin = vdmin - truncv
    tmax = vdmax + truncv
    if (tmin < min(vd)) or (tmax > max(vd)):
        check = False
    else:
        check = True
    return check
#%%
def get_abs_bound(i2,ID,vd,rsd,wtm,rad,sbh,red,rtd,truncv,vdmin,vdmax):
    '''
    This function gets the values of HITRAN lines within a provided range
    of wavenumbers.
    '''
    tmin = vdmin - truncv
    tmax = vdmax + truncv
    if (tmin < min(vd)) or (tmax > max(vd)):
        print('Min. truncation wavelength given [microns]: ',wvn2wvl(tmin))
        print('Min. lineparam wavelength [microns]: ',wvn2wvl(min(vd)))
        print('Max. truncation wavelength [microns]: ',wvn2wvl(tmax))
        print('Max. lineparam wavelength [microns]: ',wvn2wvl(max(vd)))
        raise ValueError('truncation exceeds available line parameters extent.')

    i2 = i2[(vd>tmin)&(vd<tmax)]
    ID = ID[(vd>tmin)&(vd<tmax)]       
    rsd = rsd[(vd>tmin)&(vd<tmax)]
    wtm = wtm[(vd>tmin)&(vd<tmax)]
    rad = rad[(vd>tmin)&(vd<tmax)]
    sbh = sbh[(vd>tmin)&(vd<tmax)]
    red = red[(vd>tmin)&(vd<tmax)]
    rtd = rtd[(vd>tmin)&(vd<tmax)]
    vd = vd[(vd>tmin)&(vd<tmax)]
        
    ntot = len(i2) # this variable passes the number of lines in the requested range
    
    return i2,ID,vd,rsd,wtm,rad,sbh,red,rtd,ntot
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
def salpha2py(ID,vd,rsd,rad,red,rtd,ntot,pres,temp):
    '''
    This function requires the FORTRAN script 'salpha' wrapped with f2py.
    
    It takes in the line params and outputs the line shape.
    '''
    ison = np.zeros(nvMAX)
    vdat = np.zeros(nvMAX)
    sdat = np.zeros(nvMAX)
    adat = np.zeros(nvMAX)
    edat = np.zeros(nvMAX)
    tcof = np.zeros(nvMAX) 
    
    ison[:ntot] = ID
    vdat[:ntot] = vd
    sdat[:ntot] = rsd
    adat[:ntot] = rad
    edat[:ntot] = red
    tcof[:ntot] = rtd 
     
    ison = np.reshape(ison,nvMAX,order='F')
    vdat = np.reshape(vdat,nvMAX,order='F')
    sdat = np.reshape(sdat,nvMAX,order='F')
    adat = np.reshape(adat,nvMAX,order='F')
    edat = np.reshape(edat,nvMAX,order='F')
    tcof = np.reshape(tcof,nvMAX,order='F')
    
    # reshape arrays
    sfin,alphal,alphad1 = salpha.salpha(istout,ison,vdat,sdat,adat,edat,
                                        tcof,ntot,pres,temp,ngas2imol)
    
    return sfin,alphal,alphad1
#%%
def spectrum2py(ID,vd,rsd,rad,red,rtd,ntot,nnv,pres,temp,specv):
    ''' 
    This function requires the FORTRAN script 'spectrum' wrapped with f2py.
    
    Calculating the absorption coefficient spectrum 
        using the line intensities "sfin", the Lorentz 
        halfwidth "alphaL", and the Doppler halfwidth
        (without the wavenumber) "alphaD1".
    
     spec(nvMAX) : the k-values at all gridpoints inside [vmin,vmax]
     Original: FORTRAN version written by Daphne Stam, TU Delft.
     Python: Written by Gourav Mahapatra, TU Delft
    '''
    # call the salpha function to get sfin,alphal,alphad1 variables
    sfin,alphal,alphad1 = salpha2py(ID,vd,rsd,rad,red,rtd,ntot,pres,temp)
    # rearrange the arrays that are not in FORTRAN format. 
    vdat = np.zeros(nvMAX)
    specv1 = np.zeros(nvMAX) # this is beaically the fortran format specv.
    
    vdat[:ntot] = vd
    specv1[:nnv] = specv
    
    vdat = np.reshape(vdat,nvMAX,order='F')
    specv1 = np.reshape(specv1,nvMAX,order='F')
     
    out,spec = spectrum.spectrum(istout,ntot,vdat,sfin,alphal,alphad1,specv1,nnv)
 
    return spec[:nnv]
#%%
def gaussian_IRF(x,sig=0.5,mu=0):
    ''' 
    In ckdis.py:
        
    Given a list of x-values, calculate Gaussian dist. 
    "sig" can be used to specify the std.dev. of the function in nm.
    '''
    sigv = abs(wvl2wvn(1.000)-wvl2wvn(1.000-sig*1e-3)) 
    p_x = (1/np.sqrt(2*np.pi*sigv**2))*np.exp(-((x-mu)**2/(2*sigv**2)))
    return p_x
#%%
def ksort(spec,specv,waven0,sig,irf_type='gaussian',useirf='False'):
    '''
    This function takes as input "spec" or absorption spectra 
    corresponding to a wavenumber interval "specv" and outputs
    k-distributed absorption values. 
    
    irf_type: 'spicav','gaussian'
    
    Outputs:
        abs_spec_sorted: sorted absorption spectra
           
    '''
    args_for_sort = np.argsort(spec)
    
    # sort the array
    abs_spec_sorted = np.zeros(len(spec)+1)
    abs_spec_sorted[1:] = spec[args_for_sort]
    
    # get the xwidths and x-scale
    xwidths = np.ones(len(specv))* abs(specv[0]-specv[1])
    xscale = []
    x = 0
    xscale.append(x)
    for i in range(len(xwidths)):
        x = x + xwidths[i]
        xscale.append(x)
    xscale = xscale/sum(xwidths)
    
    # get the IRF corresponding to the given wavenumber interval
    if irf_type == 'spicav': 
        irf_intp,_,dnu,_,_ = spicav_irf(0.0)
        specv_irf = specv - waven0
        irf = irf_intp(specv_irf)
        #plt.plot(specv_irf,irf_intp(specv_irf))    # this plot shows the SPICAV-IRF for a chosen sigma
    if irf_type == 'gaussian': 
        irf = gaussian_IRF(specv,sig=sig,mu=waven0)
    
    # change the widths as per the IRF after scaling the IRF.
    xwidths_conv = xwidths*(irf/irf.max())
    
    xwidths_conv_sorted = xwidths_conv[args_for_sort]
    
    # calculate the new xscale for the xwidths_conv
    xscale_conv_sorted = []
    x2 = 0
    xscale_conv_sorted.append(x2)
    for i in range(len(xwidths_conv_sorted)):
        x2 = x2 + xwidths_conv_sorted[i]
        xscale_conv_sorted.append(x2)
    xscale_conv_sorted = xscale_conv_sorted/sum(xwidths_conv_sorted)
    
#    plt.figure()
#    plt.semilogy(xscale,abs_spec_sorted,label='original spectra')
#    #plt.semilogy(xscale_conv,abs_spec_sorted,label='response applied to sorted spectra')
#    plt.semilogy(xscale_conv_sorted,abs_spec_sorted,label='response applied')
    
    return xscale,xscale_conv_sorted,abs_spec_sorted,irf_type
    
#%%    
def gauleg(nump):
    '''
    Get the Gauss points and weights by selecting the number of points (np)
    for Gauss-Legendre integration.
    
    Here range of gp is changed to [0,1] !!!
    '''
    x,w = np.polynomial.legendre.leggauss(nump)
    
    # set the variables for rescaling the interval into [a,b]
    a = 0
    b = 1
    gp = 0.5*(x + 1)*(b - a) + a
    
    gw = w
    
    # rescale the weights
    print("gauleg: The final sum after integration has to be scaled with (b-a)/2!!!")
    print("where [a,b]=[0,1].")
    
    return gp,gw
    
#%%
def kspec(wav,pres,temp,gp,sig=1.0,useirf=False,irf_type='gaussian'):
    '''
    Make k-dis spectra over a range of central wavelengths.
    '''
    # define a list of central wavelengths (in microns)
    #wav0_list = np.arange(wavmin,wavmax,step) 
    wav0_list = wav    
        
    # define the array to store k-distributed spectra per layer
    layer_kdis = []
    
    # load the whole absorption lines file
    if line_source == 'local':
        i2_a,ID_a,vd_a,rsd_a,wtm_a,rad_a,sbh_a,red_a,rtd_a,ntot_a = rdabsfile()
    else:
        i2_a,ID_a,vd_a,rsd_a,wtm_a,rad_a,sbh_a,red_a,rtd_a,ntot_a = rdHAPIparam()
    
    for wavel0 in wav0_list: 
        waven0 = wvl2wvn(wavel0) # in cm-1
        #print("Wavelength = , {:.5f} ".format(wavel0))
        
        # get the slitfunction parameters for each wav0
        nnv,truncv,specv,speci,vmin,vmax,tmin,tmax = slitfunction(wavel0,sig,nw,truncw)
        
        # check if the wavelength has absorption lines within range
        check = abs_lim_check(vmin,vmax,truncv,vd_a)
        
        if check:
            # get the values of line parameters within range [vmin-truncv,vmax+truncv] 
            i2,ID,vd,rsd,wtm,rad,sbh,red,rtd,ntot= \
                    get_abs_bound(i2_a,ID_a,vd_a,rsd_a,wtm_a,rad_a,sbh_a,red_a,rtd_a,truncv,vmin,vmax)
            # get the abs. spectra for a given p,T
            spec = spectrum2py(ID,vd,rsd,rad,red,rtd,ntot,nnv,pres,temp,specv)
        else:
            # pass empty spec to ksort.
            spec = np.zeros(specv.shape)
        
#        # get the k-sorted spectra
        x,x_irf,abs_spec,irf_type = ksort(spec,specv,waven0,sig,irf_type=irf_type,useirf=useirf)
#        # make an interpolation box
        if useirf == False:
            x_intp = interpolate.interp1d(x,abs_spec)
        else:
            x_intp = interpolate.interp1d(x_irf,abs_spec)
            print('Instrument response is being used!',' IRF type: ',irf_type)
#        # find the values of kdis at given Gauss points
        abs_gp = x_intp(gp)
#        # this should have a total len()=no.ofwavel0
        layer_kdis.append(abs_gp)
        
    return np.array(layer_kdis)    
#%%    
def sigAbs2bmAbs(p_list,t_list,atm_kdis,wav,write=False):
    '''
    This function converts sig_abs into bmabs for Venus atmosphere.
    
    It actually identifies the absorption cross-section spectra 
    at the given atmospheric boundaries and finds an average.
    
    Input:
        pressure (bars) and temperature (K)
    '''
    # add zero elements to the p_list and t_list
    p_list = list(p_list)
    t_list = list(t_list)
#    p_list.append(0.0)
#    t_list.append(0.0)
    
    if len(p_list) != len(atm_kdis):
        raise ValueError('sigAbs2bmAbs: kdis array and atmosphere array mismatch!')
    
    # convert back into array
    p_list = np.array(p_list)
    t_list = np.array(t_list)
    
    # Mean molecular weight (kg/mol)
    m_gas = 44.01e-3/6.023e23  #CO2 kg/molecule
    g_ven = 8.87               #m/s2
    
    Nd_list = []
    bmabs_kdis = []
    for i in range(len(p_list)-1):
        p_bottom = p_list[i]*1.e5
        p_top = p_list[i+1]*1.e5
        
        # find the k-distributed sig_abs of bottom and top layers
        layer_kdis_bottom = atm_kdis[i,:,:]
        layer_kdis_top = atm_kdis[i+1,:,:]
        layer_avg = np.average([layer_kdis_bottom,layer_kdis_top],axis=0)
        
        if p_bottom < p_top:
            raise ValueError('Bottom pressure lower than top pressure!')
            
        # Compute the Ngas
        Nd = (p_bottom - p_top)/(m_gas*g_ven) #/m2
        Nd_list.append(Nd)
        
        # compute bmabs for the layers average
        bmabs = layer_avg*Nd*1e-4
        bmabs_kdis.append(bmabs)
    
    bmabs_kdis = np.array(bmabs_kdis)
     
    if write == True:    
        """
        Write the bmabs values into a file for DAP code input
        """
        np.save('/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/fullKDIS/bmabs_kd.npy',bmabs_kdis)
        babsdat = open('/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/fullKDIS/bmabs.out','w')
        babsdat.write('Number of wavelength values:  \n')
        babsdat.write(str(wav.shape[0])+' \n')
        babsdat.write('Number of layers:  \n')
        babsdat.write(str(bmabs_kdis.shape[0])+' \n')
        babsdat.write('Number of gauss points:  \n')
        babsdat.write(str(bmabs_kdis.shape[2])+' \n')
        babsdat.write('lamb_no    gausspoint    layer    bmabs_kd   lambda (microns)\n')
        # file writing format:
        # wavelength_no.  gausspoint  layer    bmabs_kd
        for i in range(bmabs_kdis.shape[1]): # wavelength loop
            for j in range(bmabs_kdis.shape[2]): # gauss point loop
                for k in range(bmabs_kdis.shape[0]): # layer loop
                    babsdat.write('   '+str('{:2.0f}').format(i+1))
                    babsdat.write('    '+str('{:2.0f}').format(j+1))
                    babsdat.write('    '+str('{:2.0f}').format(k+1))
                    babsdat.write('    '+str('{:12.10E}').format(bmabs_kdis[k,i,j]))
                    babsdat.write('   '+str('{:8.4f}').format(wav[i])+'\n')
        babsdat.close()
    
    return np.array(bmabs_kdis)
    

#%%
def sigAbs2bmAbs_lbl(p_list,t_list,atm_kdis,wav,write=False):
    '''
    This function converts sig_abs into bmabs for Venus atmosphere,
    adapted for line-by-line calculation.
    
    It actually identifies the absorption cross-section spectra 
    at the given atmospheric boundaries and finds an average.
    
    Input:
        pressure (bars) and temperature (K)
    '''
    # add zero elements to the p_list and t_list
    p_list = list(p_list)
    t_list = list(t_list)
#    p_list.append(0.0)
#    t_list.append(0.0)
    
    if len(p_list) != len(atm_kdis):
        raise ValueError('sigAbs2bmAbs: kdis array and atmosphere array mismatch!')
    
    # convert back into array
    p_list = np.array(p_list)
    t_list = np.array(t_list)
    
    # Mean molecular weight (kg/mol)
    m_gas = 44.01e-3/6.023e23  #CO2 kg/molecule
    g_ven = 8.87               #m/s2
    
    Nd_list = []
    bmabs_kdis = []
    for i in range(len(p_list)-1):
        p_bottom = p_list[i]*1.e5
        p_top = p_list[i+1]*1.e5
        
        # find the k-distributed sig_abs of bottom and top layers
        layer_kdis_bottom = atm_kdis[i,:,]
        layer_kdis_top = atm_kdis[i+1,:,]
        layer_avg = np.average([layer_kdis_bottom,layer_kdis_top],axis=0)
        
        if p_bottom < p_top:
            raise ValueError('Bottom pressure lower than top pressure!')
            
        # Compute the Ngas
        Nd = (p_bottom - p_top)/(m_gas*g_ven) #/m2
        Nd_list.append(Nd)
        
        # compute bmabs for the layers average
        bmabs = layer_avg*Nd*1e-4
        bmabs_kdis.append(bmabs)
    
    bmabs_kdis = np.array(bmabs_kdis)
     
    if write == True:    
        """
        Write the bmabs values into a file for DAP code input
        """
        np.save('/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/LineByLine0/bmabs_lbl.npy',bmabs_kdis)
        np.save('/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/LineByLine0/wav.npy',wav)
        babsdat = open('/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/LineByLine0/bmabs.out','w')
        babsdat.write('Number of wavelength values:  \n')
        babsdat.write(str(wav.shape[0])+' \n')
        babsdat.write('Number of layers:  \n')
        babsdat.write(str(bmabs_kdis.shape[0])+' \n')
        #babsdat.write('Number of gauss points:  \n')
        #babsdat.write(str(bmabs_kdis.shape[2])+' \n')
        babsdat.write('lamb_no    layer    bmabs   lambda (microns)\n')
        # file writing format:
        # wavelength_no.  gausspoint  layer    bmabs_kd
        for i in range(bmabs_kdis.shape[1]): # wavelength loop
                for k in range(bmabs_kdis.shape[0]): # layer loop
                    babsdat.write('   '+str('{:2.0f}').format(i+1))
                    #babsdat.write('    '+str('{:2.0f}').format(j+1))
                    babsdat.write('    '+str('{:2.0f}').format(k+1))
                    babsdat.write('    '+str('{:12.10E}').format(bmabs_kdis[k,i]))
                    babsdat.write('   '+str('{:8.4f}').format(wav[i])+'\n')
        babsdat.close()
    
    return np.array(bmabs_kdis)
    
    
    
    
#%%
'''
Main/testing area..............................................................
'''
# the window around wavel0 (in um!!)
global sigma
sigma = 100.0  * 1.e-3      # this is for general purpose (should be 1 nm for the absorption paper simulations!!)
#sigma = 30 * 1.e-3         # this is for SPICAV-IRF # 30 nm in microns!!! SPICAV-IRF had a broad response.

# desired resolution of spectra (in points per um!!)
global nw
nw = 500000     # this is the preferred resolution!

# the number of integration points across kd spectrum
nump = 20

# get the gp,gw values for future integration
gp,gw = gauleg(nump)

# truncation window (in um!!)
# this variable determines the amount of extra wavelength to be considered at the edges 
# of wavmin and wavmax.
global truncw
truncw = 0.2
#%%
'''this is where we decide whether to test or compute a whole bmabs!'''
tester = False
recompute = True
linebyline = False

#%%  
''' This is the old 38 layer model used by Ignatiev et al.'''
# initialize the pressure and temperature arrays (this atm_profile file is different)
atm_profile2 = np.loadtxt('/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/'\
                          +'vertical_profile_ignatiev_38lays.htp',skiprows=2)
p = atm_profile2[:,2]/1.e3
t = atm_profile2[:,1]
alt = atm_profile2[:,0]

# make the interpolation objects
int_temp = interpolate.interp1d(alt,t)
int_pres = interpolate.interp1d(alt,p)

#%% 
'''
This is where the 71 layer atmosphere is created!!!
'''
# get the new atmosphere with 1 km resolution above 40 km.
alt_low = np.array([     0., 4., 8., 12., 16., 20., 
           24., 28., 32., 36., 40.])
alt_high = np.arange(41,101,1)     
alt = np.zeros(len(alt_low)+len(alt_high))
alt[0:len(alt_low)] = alt_low
alt[len(alt_low):] = alt_high

# get the nwe p, t as per the new alt
# interpolate
t = int_temp(alt)
p = int_pres(alt)

p_list = np.zeros(len(p)+1)
t_list = np.zeros(len(t)+1)

p_list[:len(p)] = p
t_list[:len(t)] = t
t_list[-1] = t_list[-2]
#%% for testing or plotting a small interval.
if tester:
    ''' This is a test function for checking the ckdis working'''
    # define the central wavelength (tester)
    wavel0 = 1.45
    waven0 = wvl2wvn(wavel0)
    
    # get all the slitfunction parameters (tester)
    nnv,truncv,specv,speci,vmin,vmax,tmin,tmax = slitfunction(wavel0,sigma,nw,truncw)

    # read the absorption lines file by Rothman et al., (2004) (tester)
    if line_source == 'local':
        i2_a,ID_a,vd_a,rsd_a,wtm_a,rad_a,sbh_a,red_a,rtd_a,ntot = rdabsfile()    
        i2,ID,vd,rsd,wtm,rad,sbh,red,rtd,ntot= \
                        get_abs_bound(i2_a,ID_a,vd_a,rsd_a,wtm_a,rad_a,sbh_a,red_a,rtd_a,truncv,vmin,vmax)
    else:
        i2_a,ID_a,vd_a,rsd_a,wtm_a,rad_a,sbh_a,red_a,rtd_a,ntot = rdHAPIparam()
        i2,ID,vd,rsd,wtm,rad,sbh,red,rtd,ntot= \
                        get_abs_bound(i2_a,ID_a,vd_a,rsd_a,wtm_a,rad_a,sbh_a,red_a,rtd_a,truncv,vmin,vmax)                    
    print('Number of lines within the truncation limits ,', ntot)

    # get the absorption spectra for a given (set) of p,T (tester)
    lev = -1
    temp = atm_profile[lev,2]      # in K
    pres = atm_profile[lev,1]      # in bars
    
    spec = spectrum2py(ID,vd,rsd,rad,red,rtd,ntot,nnv,pres,temp,specv)
    #spec = np.zeros(specv.shape)
    
    #plt.figure(figsize=(12,8))
    plt.subplot(211)
    plt.semilogy(wvn2wvl(specv),spec,ls=ls.linestyles['densely dashed'],lw=1,label='{:.6f}'.format(pres),c='k')
    plt.grid()
    plt.xlabel('Wavelength ($\mu m$)',fontsize='large')
    plt.ylabel('$\sigma_{abs} (cm^{-2})$',fontsize='large')
    
    # call the KDIS function (tester)
    x0,x1,ksp,_ = ksort(spec,specv,waven0,0.005,irf_type='spicav')
    plt.subplot(212)
    #plt.semilogy(x1,ksp)
    plt.semilogy(x0,ksp,ls=ls.linestyles['solid'],lw=1,label='{:2.6f} bar'.format(pres),c='k')    # this plots the ck-dis without the intstrument response
    plt.xlabel('x-scale',fontsize='large')
    plt.ylabel('$\sigma_{abs} (cm^{-2})$',fontsize='large')
    plt.grid()
    plt.legend(loc=4)
    #x0,x2,ksp = ksort(spec,specv,waven0,0.1)
    #x0,x3,ksp = ksort(spec,specv,waven0,0.2)
    
    #get the gaussian response (tester)
    #irf1 = gaussian_IRF(specv,sig=0.05,mu=waven0)
    #irf2 = gaussian_IRF(specv,sig=0.1,mu=waven0)
    #irf3 = gaussian_IRF(specv,sig=0.2,mu=waven0)

    # test for the spectra output (tester)
#    spec = []
#    for i in range(len(p_list)):
#        spec.append(spectrum2py(ID,vd,rsd,rad,red,rtd,ntot,nnv,p_list[i],t_list[i],specv))
#    spec = np.array(spec)
#%% for the line-by-line case i.e. without any ck-distribution!
if linebyline:
    '''
    This section computes bmabs for only line-by-line spectra.
    '''
    #wav = np.load('/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/fullKDIS/wav.npy')
    
    ''' Form a custom geometry function for the input into fullKDIS'''
    ngeos = 91
    phase = np.linspace(0,90,ngeos)
    sza = phase
    emission = np.zeros(phase.shape)
    phi = np.zeros(phase.shape)
    beta = np.zeros(phase.shape)
    
    geos_file = open('/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/geos.in','w')
    geos_file.write('# this file contains the geometries obtained from VEx SPICAV IR dataset\n')
    geos_file.write('# number of available geometries (ngeos):\n')
    geos_file.write(str(ngeos))
    geos_file.write('#    alpha    theta0    theta    phi    beta \n')
    for i in range(ngeos):
        geos_file.write('       {:2.2f}      {:2.2f}      {:2.2f}      {:2.2f}      {:2.2f} \n'.format(phase[i],sza[i],emission[i],phi[i],beta[i]))
    geos_file.close()
    
    ''' calculate the molecular absorption optical thickness of each layer'''

    lbl_spec = []
    
    for i in range(len(p_list)):
        print("Layer ",i+1,"...")
        # start the timer
        t0 = time.time()    
        pres = p_list[i]
        temp = t_list[i]
        
        # define the central wavelength (tester)
        wavel0 = 1.45
        # define the sigma (this is special case for line-by-line)
        sigma = 100.0  * 1.e-3
        
        waven0 = wvl2wvn(wavel0)
    
        # get all the slitfunction parameters (tester)
        nnv,truncv,specv,speci,vmin,vmax,tmin,tmax = slitfunction(wavel0,sigma,nw,truncw)
        
        # convert the specv into wavelength
        wav = wvn2wvl(specv)

        # read the absorption lines file by Rothman et al., (2004) (tester)
        if line_source == 'local':
            i2_a,ID_a,vd_a,rsd_a,wtm_a,rad_a,sbh_a,red_a,rtd_a,ntot = rdabsfile()    
            i2,ID,vd,rsd,wtm,rad,sbh,red,rtd,ntot= \
                        get_abs_bound(i2_a,ID_a,vd_a,rsd_a,wtm_a,rad_a,sbh_a,red_a,rtd_a,truncv,vmin,vmax)
        else:
            i2_a,ID_a,vd_a,rsd_a,wtm_a,rad_a,sbh_a,red_a,rtd_a,ntot = rdHAPIparam()
            i2,ID,vd,rsd,wtm,rad,sbh,red,rtd,ntot= \
                        get_abs_bound(i2_a,ID_a,vd_a,rsd_a,wtm_a,rad_a,sbh_a,red_a,rtd_a,truncv,vmin,vmax)                    
        print('Number of lines within the truncation limits ,', ntot)
        
        spec = spectrum2py(ID,vd,rsd,rad,red,rtd,ntot,nnv,pres,temp,specv)
        
        lbl_spec.append(spec)
        t_tot = (time.time() - t0)*len(p_list)
        print("Time remaining provided there are ",len(p_list),' layers is: ',(t_tot-i*(time.time() - t0))/60,' mins.')
    lbl_spec = np.array(lbl_spec)
    
    # reverse the order of the storing of wav and bmabs:
    wavnew = []
    lbl_specnew = []
    for i in range(len(wav)-1,0,-1):
        wavnew.append(wav[i])
        lbl_specnew.append(lbl_spec[:,i])
        
    # convert to array
    wav = np.array(wavnew)
    lbl_spec = np.transpose(np.array(lbl_specnew))
    
    bmabs_lbl = sigAbs2bmAbs_lbl(p_list,t_list,lbl_spec,wav,write=True)
    
#%% this is where the ck-distributed spectra is created 
if not tester and recompute: 
    ''' This assumes you are creating a KDIS spectra'''   
    # initialize the wavelength list here!!
    #wavmin = 0.5
    #wavmax = 2.005
    #wavmin = 1.4   # in microns
    #wavmax = 1.5   # in microns
    #step = sigma/2
    #step = 0.005    # in microns
    wav = np.load('/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/fullKDIS/wav.npy')
    #wav = np.arange(wavmin,wavmax,step)

    # write the wav file
    #np.save('/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/fullKDIS/wav5nm_500to2000.npy',wav)
#    wavdat = open('/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/wav.dat','w')
#    for i in range(len(wav)):
#        wavdat.write(str(i+1)+'   ')
#        wavdat.write(str('%.4f   '%wav[i])+'\n')
#    wavdat.close() 
    ''' Form a custom geometry function for the input into fullKDIS'''
    ngeos = 91
    phase = np.linspace(0,90,ngeos)
    sza = phase
    emission = np.zeros(phase.shape)
    phi = np.zeros(phase.shape)
    beta = np.zeros(phase.shape)
    
    geos_file = open('/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/geos.in','w')
    geos_file.write('# this file contains the geometries obtained from VEx SPICAV IR dataset\n')
    geos_file.write('# number of available geometries (ngeos):\n')
    geos_file.write(str(ngeos))
    geos_file.write('#    alpha    theta0    theta    phi    beta \n')
    for i in range(ngeos):
        geos_file.write('       {:2.2f}      {:2.2f}      {:2.2f}      {:2.2f}      {:2.2f} \n'.format(phase[i],sza[i],emission[i],phi[i],beta[i]))
    geos_file.close()
    
    ''' calculate the molecular absorption optical thickness of each layer'''

    atm_kdis = []
    for i in range(len(p_list)):
        print("Layer ",i+1,"...")
        # start the timer
        t0 = time.time()    
        pres = p_list[i]
        temp = t_list[i]
        layer_kdis = kspec(wav,pres,temp,gp,sig=sigma,useirf=False,irf_type = 'spicav')
        atm_kdis.append(layer_kdis)
        t_tot = (time.time() - t0)*len(p_list)
        print("Time remaining provided there are ",len(p_list),' layers is: ',(t_tot-i*(time.time() - t0))/60,' mins.')
    atm_kdis = np.array(atm_kdis)
        
    bmabs_kdis = sigAbs2bmAbs(p_list,t_list,atm_kdis,wav,write=True)

