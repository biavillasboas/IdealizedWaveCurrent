import numpy as np

import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib
import xarray as xr
import cmocean.cm as cmo

import os

matplotlib.rcParams['axes.linewidth'] = 0.8
matplotlib.rcParams['axes.edgecolor']='black'
matplotlib.rcParams['font.family'] = 'serif'
# This requires a working LaTeX installation. If you don't have that,
# you may comment this line
matplotlib.rcParams['text.usetex']=True

figpath = "../figs/"
expnames = ['standard', 'double_energy']

n = 20
slope = 2.5
period = 10.3

lin_cur = np.linspace(0, 0.7, 30)
lin_hs = np.linspace(0.8,1.2, 30)
letc = ['A', 'B']
leth = ['C', 'D']
bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=1.5, alpha=.9)
seed = 15 
fig, axes = plt.subplots(2,2,sharex='all', sharey='all', figsize=(6,6))
for i, exp in enumerate(expnames):
    path = "../data/synthetic/%s/gridded/" %exp
    if exp=='standard':
        div=0.0
        std = 0.01
    else:
        div=0.5
        std = 0.02

    cur_filename = "K%sA%sT%sS%s_cur.nc" %(slope, div, period, std)
    cur_data = os.path.join(path, cur_filename)
    dsc = data = xr.open_dataset(cur_data)
    x = dsc.x*1e-3
    u = dsc['ucur'][seed]
    v = dsc['vcur'][seed]
    U = (u**2 +   v**2)**.5

    hs_filename = "K%sA%sT%sS%s_hs.nc" %(slope, div, period, std)
    hs_data = os.path.join(path, hs_filename)
    dsh = xr.open_dataset(hs_data)
    hs = dsh.hs[seed,]
    axc = axes[0][i]
    axh = axes[1][i]
    csc = axc.contourf(x, x, U, lin_cur,  extend='both', cmap='GnBu_r')
    axc.quiver(x[::n], x[::n], u[::n,::n]/U[::n,::n], v[::n,::n]/U[::n,::n], color='black', width=0.004, scale=25)
    csh = axh.contourf(x, x, hs, lin_hs, extend='both', cmap=cmo.curl)

    axc.set_xticks([50, 250, 450])
    axh.set_xticklabels([50, 250, 450])
    axh.set_yticks([50, 250, 450])
    axh.set_yticklabels([50, 250, 450])
    axc.set_yticks([50, 250, 450])
    axc.set_yticklabels([50, 250, 450])
    plt.subplots_adjust(wspace=0.0, hspace=0.0)
    axh.text(40, 550, '%s' %leth[i], fontsize=14, bbox=bbox_props)
    axc.text(40, 550, '%s' %letc[i], fontsize=14, bbox=bbox_props)

    # rect = [left, bottom, width, height]
    cax1 = plt.axes([0.92, 0.14, 0.015, 0.3])
    cax2 = plt.axes([0.92, 0.53, 0.015, 0.3])
    cbar2 = plt.colorbar(csc, cax=cax2, extend='both', ticks=[0, 0.2, 0.4, 0.6], orientation='vertical',format = '%.1f')
    cbar2.ax.xaxis.set_ticks_position('bottom')
    cbar2.ax.xaxis.set_label_position('bottom')
    cbar2.set_label('$U$ [m/s]',fontsize=14, labelpad=14)
    cbar2.ax.tick_params(labelsize=12)
    cbar2.update_ticks()

    cbar1 = plt.colorbar(csh,cax=cax1, extend='both', orientation='vertical',format = '%.1f')
    cbar1.ax.tick_params(labelsize=12)
    cbar1.ax.xaxis.set_ticks_position('top')
    cbar1.ax.xaxis.set_label_position('top')
    cbar1.set_label('$H_s$ [m]', fontsize=14,  labelpad=14)
    tick_locator = ticker.MaxNLocator(nbins=4)
    cbar1.locator = tick_locator
    cbar1.update_ticks()
fig.text(0.05, 0.5, 'distance [km]', ha='center', va='center', rotation='vertical', fontsize=14)
fig.text(0.51, 0.03, 'distance [km]', ha='center', va='center', fontsize=14)
plt.subplots_adjust(wspace=0.03, hspace=0.03)
figname= "fig08.png"
plt.savefig(figpath+figname, dpi=150, bbox_inches='tight')
