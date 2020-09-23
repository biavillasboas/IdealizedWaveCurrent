import sys 
sys.path.append('../tools')

import os
import numpy as np
import glob
import xarray as xr
import xrft
from netCDF4 import Dataset, num2date

from scipy import stats
import pandas as pd
from spectra_tools import *
from stochastic_flow_tools import *

def get_parameters(data):
    """ Extracts parameters from netcdf attributes
    """
    parameters_dict = data.attrs 
    return parameters_dict

def basic_stats(var, data, ni, nf):
    """ Computes basic statistics of gridded output
    """
    output = {}
    dx = float(data[var].x[1] - data[var].x[0])
    da = data[var][ni:nf, ni:nf]
    output['%s_var' %var] = float(da.var())
    output['%s_std' %var] = float(da.std())
    output['%s_median'%var] = float(da.median())
    output['%s_mean'%var] = float(da.mean())
    output['%s_skew'%var] = stats.skew(da.values.flatten())
    output['%s_kurt'%var] = stats.kurtosis(da.values.flatten())
    fy, fx = np.gradient(da.values, dx)
    grad2 = fx**2 + fy**2
    output['%s_grad'%var] = np.sqrt(np.mean(grad2))
    return output

def eke_spectral_slope(u, v):
    """ Computes spectral slope of the KE wavenumber spectrum between 10km and 100km
    """
    u_iso = xrft.isotropic_powerspectrum(u, dim=['x','y'], detrend='linear', window=True)
    v_iso = xrft.isotropic_powerspectrum(v, dim=['x','y'], detrend='linear', window=True)
    Ei = 0.5*(u_iso + v_iso)
    df = np.diff(Ei.freq_r[1:]).mean()
    diffusivity = np.sum(Ei.freq_r*Ei*df)
    q1=1/100e3
    q2 = 1/10e3
    ind = ((Ei.freq_r>=q1) & (Ei.freq_r<=q2)).values
    s, _, _, _, _ = stats.linregress(np.log10(Ei.freq_r[ind].values), np.log10(Ei[ind].values))
    return abs(s), float(diffusivity) 


def flow_stats(var, data, ni, nf):
    output = {}
    u = data.ucur[ni:nf, ni:nf]
    v = data.vcur[ni:nf, ni:nf]
    dx = (data.x.values[1] - data.x.values[0])
    up, vp, us, vs = helm_decompose_flow(u, v, dx)
    vort, div = flow_gradients(u, v, dx)
    eke = 0.5*np.mean((u**2 + v**2))
    ekep = 0.5*np.mean((up**2 + vp**2))
    ekes = 0.5*np.mean((us**2 + vs**2))
    output['div_ratio'] = float(ekep/eke)
    output['vorticity'] = np.sqrt(np.mean(vort**2))
    output['divergence'] = np.sqrt(np.mean(div**2))
    output['EKE'] = float(eke)
    output['EKE_phi'] = float(ekep)
    output['EKE_psi'] = float(ekes)
    output['EKE_slope'], __ = eke_spectral_slope(u, v) 
    output['EKE_psi_slope'], output['diffusivity'] = eke_spectral_slope(us, vs) 
    output['EKE_phi_slope'],  __ = eke_spectral_slope(up, vp) 
    return output

def hs_spectral_slope(var, data, ni, nf):
    """ Computes spectral slope of the Hs wavenumber spectrum between 10km and 100km
    """
    output = {}
    hs = data[var][ni:nf, ni:nf]
    hs_hat2 = xrft.power_spectrum(hs, dim=['x','y'], detrend='linear', window=True)
    nk = hs_hat2.shape[-1]
    hsy = hs_hat2.sum(dim='freq_x')[nk//2:]*hs_hat2.freq_x_spacing
    hsx = hs_hat2.sum(dim='freq_y')[nk//2:]*hs_hat2.freq_y_spacing
    q1=1/100e3
    q2 = 1/10e3
    ind = ((hsy.freq_y>=q1) & (hsy.freq_y<=q2)).values
    sy, intercepty, r_valuey, p_valuey, std_erry = stats.linregress(
            np.log10(hsy.freq_y[ind].values), np.log10(hsy[ind].values))
    sx, _, _, _, _ = stats.linregress(np.log10(hsx.freq_x[ind].values), np.log10(hsx[ind].values))
    output['hs_yslope'] = abs(sy) 
    output['hs_xslope'] = abs(sx) 
    return output


def analize_member(filename, var, diagnostic_functions):
    """ Main function to process the gridded data from synthetic runs
    """
    
    data = xr.open_dataset(filename)
    ni = np.where(data.x==200000)[0][0]
    nf = np.where(data.x==600000)[0][0]
    parameters_dict = get_parameters(data)
    output = []
    for seed in data.seed:
        member = data.sel(seed=seed)
        output.append({})
        for dfunc in diagnostic_functions:
            tmp = dfunc(var, member, ni, nf)
            for pm in tmp:
                output[-1][pm] = tmp[pm]
    ds = xr.Dataset()
    ds.coords['alpha'] = ('alpha', np.atleast_1d(parameters_dict['alpha']))
    ds.coords['period'] = ('period', np.atleast_1d(parameters_dict['period']))
    ds.coords['slope'] = ('slope', np.atleast_1d(parameters_dict['slope']))
    ds.coords['seed'] = ('seed', data.seed)
    df = pd.DataFrame(output)
    for key in df.keys():
        ds[key] = xr.DataArray([[[df[key]]]], coords=(ds['alpha'], ds['slope'], ds['period'], ds['seed']))
    return ds


def analize_llc_member(filename, var, diagnostic_functions):
    """ Main function to process the gridded data from llc runs
    """
    
    data = xr.open_dataset(filename)
    ni = np.where(data.x==200000)[0][0]
    nf = np.where(data.x==600000)[0][0]
    nc = Dataset(filename, 'r')
    output = []
    for t in data.time:
        member = data.sel(time=t, drop=True)
        output.append({})
        for dfunc in diagnostic_functions:
            tmp = dfunc(var, member, ni, nf)
            for pm in tmp:
                output[-1][pm] = tmp[pm]
    df = pd.DataFrame(output)
    df['time'] = num2date(nc.variables['time'][:], nc.variables['time'].units)
    ds = df.set_index('time').to_xarray()
    return ds
