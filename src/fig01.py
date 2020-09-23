import sys
sys.path.append('../tools')

import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from matplotlib import ticker, rc
import matplotlib

import xrft
import xarray as xr

from spectra_tools import *
from stochastic_flow_tools import *

matplotlib.rcParams['xtick.labelsize'] = 12
matplotlib.rcParams['ytick.labelsize'] = 12
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['text.usetex']=True


figpath = '../figs/'

L = 500e3
dx = 2.5e3
slopes = [1.66, 3.0]
bbox_props = dict(boxstyle="round,pad=0.3", fc="w", ec="k", lw=1, alpha=.95)
st = {'1.66':'5/3','3.0':'3.0'}

lin_cur = np.linspace(0, 0.4, 30)
n=10

fig, axes = plt.subplots(2, 3, sharex='col', figsize=(12,6))
axes_s = axes[:,0]
axes_rot = axes[:,1]
axes_div = axes[:,2]
for i, slope in enumerate(slopes):

    u_psi, v_psi, u_phi, v_phi = potential_solenoidal(L=L, dx=dx, slope=slope,
                                                      kmin=2*pi/300e3, kmax=2*pi/5e3, seed=10)
    u_psi, v_psi, u_phi, v_phi, u_mix, v_mix = synthetic_flow(
                                              u_psi, v_psi, u_phi, v_phi,
                                              mean_eke=0.01, alpha=0.5)
    x = np.arange(0, L, dx)*1e-3
    y = np.arange(0, L, dx)*1e-3

    u_psi = xr.DataArray(u_psi,
            dims=('y', 'x'), coords={'y':y, 'x':x})
    v_psi = xr.DataArray(v_psi,
            dims=('y', 'x'), coords={'y':y, 'x':x})
    u_phi = xr.DataArray(u_phi,
            dims=('y', 'x'), coords={'y':y, 'x':x})
    v_phi = xr.DataArray(v_phi,
            dims=('y', 'x'), coords={'y':y, 'x':x})
    u_mix = xr.DataArray(u_mix,
            dims=('y', 'x'), coords={'y':y, 'x':x})
    v_mix = xr.DataArray(v_mix,
            dims=('y', 'x'), coords={'y':y, 'x':x})
    Us = (u_psi**2 + v_psi**2)**.5
    Up = (u_phi**2 + v_phi**2)**.5
    Um = (u_mix**2 + v_mix**2)**.5

    um_iso = xrft.isotropic_powerspectrum(u_mix, dim=['x','y'], detrend='linear', window=True)
    vm_iso = xrft.isotropic_powerspectrum(v_mix, dim=['x','y'], detrend='linear', window=True)
    Ei_mix = 0.5*(um_iso + vm_iso)

    axpsi = axes_rot[i]
    axs = axes_s[i]
    axphi = axes_div[i]

    axs.loglog(Ei_mix.freq_r, Ei_mix, color='k', label='q$^{-%s}$'%st['%s'%slope])
    axs.legend(loc='upper right', fontsize=12)
    axs.set_ylim([1e-7, 1e-1])
    axs.grid(which='both', ls='dotted')
    axpsi.contourf(x, x, Us, lin_cur, extend='both', cmap='GnBu_r')
    axpsi.quiver(x[::n], x[::n], u_psi[::n,::n]/Us[::n,::n], v_psi[::n,::n]/Us[::n,::n],
            color='black', width=0.004, scale=25)
    
    cs = axphi.contourf(x, x, Up, lin_cur, extend='both', cmap='GnBu_r')
    axphi.quiver(x[::n], x[::n], u_phi[::n,::n]/Up[::n,::n], v_phi[::n,::n]/Up[::n,::n],
            color='black', width=0.004, scale=25)
    axpsi.yaxis.set_ticklabels('')
    axphi.yaxis.tick_right()
    axphi.yaxis.set_label_position('right')
    if i == 0:
        axs.xaxis.set_ticklabels('')
cax1 = fig.add_axes([0.5, 0.9, 0.3, 0.02])
cbar1 = plt.colorbar(cs,cax=cax1, extend='both', orientation='horizontal',format = '%.1f')
cbar1.ax.xaxis.set_ticks_position('top')
cbar1.ax.xaxis.set_label_position('top')
tick_locator = ticker.MaxNLocator(nbins=4)
cbar1.locator = tick_locator
cbar1.update_ticks()
cbar1.set_label('$U$ [m/s]', fontsize=14, labelpad=14)
fig.text(0.65, 0.03, 'distance [km]', ha='center', va='center', fontsize=14)
fig.text(0.05, 0.5, r"$\tilde E(q)$ [m$^2$/s$^2$/cycles/km]", rotation='vertical',ha='center', va='center', fontsize=14)
fig.text(0.24, 0.03, 'wavenumber [cycles/km]', ha='center', va='center', fontsize=14)
fig.text(0.96, 0.52, 'distance [km]', rotation=270,ha='center', va='center', fontsize=14)
fig.subplots_adjust(wspace=0.1, hspace=0.1)
fig.savefig(figpath + 'fig01.png', dpi=150, bbox_inches='tight')
