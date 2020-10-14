import sys
sys.path.append('../tools')

import numpy as np

import matplotlib.pyplot as plt
import matplotlib
import xrft
import xarray as xr

from spectra_tools import spec_error

import os

matplotlib.rcParams['axes.linewidth'] = 0.8
matplotlib.rcParams['axes.edgecolor']='black'
matplotlib.rcParams['font.family'] = 'serif'
# This requires a working LaTeX installation. If you don't have that,
# you may comment this line
matplotlib.rcParams['text.usetex']=True

k1 = 1/200
k2 = 1/15

figpath = "../figs/"
expnames = ['standard', 'double_energy']

slopes = [1.66, 2.0, 2.5, 3.0]
period = 10.3
st = {'1.66':'5/3', '2.0':'2.0', '2.5':'2.5', '3.0':'3.0'}

bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=1, alpha=1)
i=1
plt.figure(figsize=(12, 6))
for slope in slopes:
    ax = plt.subplot(1,4,i)
    for exp in expnames:
        path = "../data/synthetic/%s/gridded/" %exp

        if exp=='standard':
            div=0.0
            std = 0.01
            ls='-'
            utype = r'$\tilde E^\psi$'
            hstype = r'$H_s^\psi$'
        else:
            div=0.5
            std = 0.02
            ls='--'
            utype = r'$\tilde E^\psi + \tilde E^\phi$'
            hstype = r'$H_s^\psi$ + $H_s^\phi$'

        cur_filename = "K%sA%sT%sS%s_cur.nc" %(slope, div, period, std)
        cur_data = os.path.join(path, cur_filename)
        dsc = data = xr.open_dataset(cur_data)
        ni = np.where(dsc.x==200000)[0][0]
        nf = np.where(dsc.x==600000)[0][0]
        u_mix = dsc['ucur'][:, ni:nf, ni:nf]
        v_mix = dsc['vcur'][:, ni:nf, ni:nf]

        um_hat2 = xrft.power_spectrum(u_mix, dim=['x','y'], detrend='linear', window=True)
        vm_hat2 = xrft.power_spectrum(v_mix, dim=['x','y'], detrend='linear', window=True)
        nk = um_hat2.shape[-1]
        ekem_hat2 = 0.5*(um_hat2 + vm_hat2)
        ekem_y = ekem_hat2.sum(dim='freq_x')[:,nk//2:]*ekem_hat2.freq_x_spacing
        N = ekem_y.seed.size
        ky = ekem_y.freq_y.values*1e3
        Em = ekem_y.mean(dim='seed').values*1e-3
        Em_l, Em_u = spec_error(Em, N, ci=0.95)
        
        hs_filename = "K%sA%sT%sS%s_hs.nc" %(slope, div, period, std)
        hs_data = os.path.join(path, hs_filename)
        dsh = xr.open_dataset(hs_data)
        hs = dsh.hs[:, ni:nf, ni:nf]

        hs_hat2 = xrft.power_spectrum(hs, dim=['x','y'], detrend='linear', window=True)
        nk = hs_hat2.shape[-1]
        hsy = hs_hat2.sum(dim='freq_x')[:,nk//2:]* hs_hat2.freq_x_spacing
        E_hsy = hsy.mean(dim='seed').values*1e-3

        ax.loglog(ky, Em, color='k', ls=ls, label='%s' %utype)
        ax.loglog(ky, E_hsy/10, color='mediumseagreen', ls=ls, label='%s' %hstype)
        ax.legend(loc='lower left', fontsize=13, framealpha=1)
        ax.axvspan(k1, k2, alpha=0.1, color='gray')
    ax.plot(np.ones(5)*3.5e-3, np.linspace(Em_l[0]*1e-3, Em_u[0]*1e-3, 5), color='gray', ls='-', linewidth=2.0)

    if i !=1:
        ax.set_yticklabels([], minor=True)
        ax.yaxis.set_ticks_position('none')
        plt.setp(ax.get_yticklabels(), visible=False)
    ax.grid('True', which='both', ls='dotted', lw=.5)
    ax.annotate('$q^{-%s}$' %st['%s'%slope], xy=(0.7,0.8), fontsize=13, xycoords="axes fraction", bbox=bbox_props)
    i+=1
    ax.set_ylim([1e-7, 1])
    ax.set_xlim([2e-3, 1e-1])
plt.annotate('Wavenumber [cycles/km]', xy=(0.43,0.03), xycoords="figure fraction", fontsize=14)
plt.annotate('PSD [m$^2$/cycle/km] or [m$^2$/s$^2$/cycle/km]', xy=(0.02,0.25), xycoords="figure fraction", fontsize=14, rotation =90)
plt.subplots_adjust(wspace=0.05, hspace=0)
figname= "fig09.png"
plt.savefig(figpath+figname, dpi=150, bbox_inches='tight')
