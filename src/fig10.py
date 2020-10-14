import sys
sys.path.append('../tools')

import glob
import numpy as np
from netCDF4 import Dataset

import cmocean.cm as cmo
import matplotlib.pyplot as plt
from matplotlib import ticker, cm
import matplotlib

import seawater as sw

from stochastic_flow_tools import flow_gradients

matplotlib.rcParams['axes.linewidth'] = 0.8
matplotlib.rcParams['axes.edgecolor']='black'
matplotlib.rcParams['font.family'] = 'serif'
# This requires a working LaTeX installation. If you don't have that,
# you may comment this line
matplotlib.rcParams['text.usetex']=True

figpath = "../figs/"
f = sw.f(31.5)

bbox_props = dict(fc="w", ec="k", lw=1, alpha=.9)
lin_hs = np.linspace(0.8,1.2,30)
lin_zeta = np.linspace(-0.92,0.91,30)
months = {'02':'Winter', '08':'Summer'} 

fig1, (axes1, axes2) = plt.subplots(2,2,sharex='all', sharey='all', figsize=(6,6))
for idx, month in enumerate(months.keys()):

    hs_file = '../data/llc4320/gridded/llc_T10.3_2012%s_hs.nc' %month
    nch = Dataset(hs_file, 'r')
    hs = nch.variables['hs'][0]
    x = nch.variables['x'][:]*1e-3
    y = nch.variables['y'][:]*1e-3

    cur_file = '../data/llc4320/gridded/llc_T10.3_2012%s_cur.nc' %month
    ncc = Dataset(cur_file, 'r')
    u = ncc.variables['ucur'][0]
    v = ncc.variables['vcur'][0]
    vort, div = flow_gradients(u, v, 2.5e3)
    zeta = vort/f

    ax1 = axes2[idx]
    ax2 = axes1[idx]
    cs1 = ax1.contourf(x, y, hs, lin_hs, cmap='cmo.curl', extend='both')
    cs2 = ax2.contourf(x, y, zeta, lin_zeta, cmap=cmo.balance, extend='both')
    ax2.text(40, 550, '%s' %months[month], fontsize=14, bbox=bbox_props)
    ax2.set_xticks([50, 250, 450])
    ax1.set_xticklabels([50, 250, 450])
    ax1.set_yticks([50, 250, 450])
    ax1.set_yticklabels([50, 250, 450])
    ax2.set_yticks([50, 250, 450])
    ax2.set_yticklabels([50, 250, 450])
    plt.subplots_adjust(wspace=0.0, hspace=0.0)
    
    # rect = [left, bottom, width, height]
    cax1 = plt.axes([0.92, 0.14, 0.015, 0.3])
    cax2 = plt.axes([0.92, 0.53, 0.015, 0.3])
    cbar2 = plt.colorbar(cs2,cax=cax2, extend='both', ticks=[-0.9, -0.45, 0, 0.45, 0.9], orientation='vertical',format = '%.2f')
    cbar2.set_ticklabels(['-0.9', '-0.45', '0', '0.45', '0.9'])
    cbar2.ax.xaxis.set_ticks_position('bottom')
    cbar2.ax.xaxis.set_label_position('bottom')
    cbar2.set_label(r'$\zeta/ f$  [ ]',fontsize=14)
    cbar2.ax.tick_params(labelsize=12)
    cbar2.update_ticks()

    cbar1 = plt.colorbar(cs1,cax=cax1, extend='both', orientation='vertical',format = '%.1f')
    cbar1.ax.tick_params(labelsize=12)
    cbar1.ax.xaxis.set_ticks_position('top')
    cbar1.ax.xaxis.set_label_position('top')
    cbar1.set_label(r'$H_s$  [m]', fontsize=14)
    tick_locator = ticker.MaxNLocator(nbins=4)
    cbar1.locator = tick_locator
    cbar1.update_ticks()
fig1.text(0.01, 0.5, 'distance [km]', ha='center', va='center', rotation='vertical', fontsize=14)
fig1.text(0.51, 0.03, 'distance [km]', ha='center', va='center', fontsize=14)
plt.subplots_adjust(wspace=0.03, hspace=0.03)
figname= "fig10.png" 
plt.savefig(figpath+figname, dpi=150, bbox_inches='tight')
