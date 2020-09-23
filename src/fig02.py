import numpy as np
from netCDF4 import Dataset

import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib
import cmocean.cm as cmo

import os

matplotlib.rcParams['axes.linewidth'] = 0.8
matplotlib.rcParams['axes.edgecolor']='black'
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['text.usetex']=True

path = "../data/synthetic/standard/gridded/"
figpath = "../figs/"

period = 10.3
slopes = [1.66, 2.0, 2.5, 3.0]
divs = [0.0, 0.4, 0.8, 1.0]
################################################
cur_lin = np.linspace(0, .5, 30)
bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=1, alpha=.9)
n = 20
nd = 250
std = 0.01
st = {'1.66':'5/3', '2.0':'2.0', '2.5':'2.5', '3.0':'3.0'}


seed_idx = 13
i=1
fig = plt.figure(figsize=(9.9,10))
for slope in slopes:
    for div in divs:
        filename = "K%sA%sT%sS%s_cur.nc" %(slope, div, period, std)
        data = os.path.join(path, filename)
        nc = Dataset(data, 'r')

        u = nc.variables['ucur'][seed_idx,:nd,:nd]
        v = nc.variables['vcur'][seed_idx, :nd, :nd]
        U = (u**2 + v**2)**.5
        x = (nc.variables['x'][:nd] - nc.variables['x'][0])*1e-3
        y = (nc.variables['y'][:nd] - nc.variables['y'][0])*1e-3

        ax = plt.subplot(4, 4,i) 
        cs = ax.contourf(x, y, U, cur_lin, extend='both', cmap='GnBu_r')
        ax.quiver(x[::n], y[::n], u[::n,::n]/U[::n,::n], v[::n,::n]/U[::n,::n], color='black', width=0.006, scale=20)

        ax.set_xticks([50, 250, 450])
        ax.set_yticks([50, 250, 450])
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top') 
        ax.tick_params(labelbottom='off',labeltop='on')

        if i not in [1, 2, 3, 4]:
            ax.set_xticklabels([])
            plt.setp(ax.get_xticklabels(), visible=False)
            plt.setp(ax.get_xticklabels(), visible=False)

        if i not in [1,5,9,13]:
            ax.set_yticklabels([])
            plt.setp(ax.get_yticklabels(), visible=False)

        ax.set_aspect('equal')
        ax.text(50, 50, '$\\alpha = %s$ \n$q^{-%s}$' %(div, st['%s'%slope]), 
                fontsize=12, bbox=bbox_props)
        i+=1
fig.text(0.05, 0.5, 'distance [km]', ha='center', va='center', rotation='vertical', fontsize=14)
fig.text(0.51, 0.95, 'distance [km]', ha='center', va='center', fontsize=14)
seed = nc.variables['seed'][seed_idx]
plt.subplots_adjust(wspace=0.0, hspace=0.0)
# rect = [left, bottom, width, height]
cax = plt.axes([0.26, 0.06, 0.5, 0.02])
ticks = [0,0.1,0.2,0.3,0.4,0.5]
cbar = plt.colorbar(cs,cax=cax, extend='both', ticks=ticks, orientation='horizontal',format = '%.1f')
cbar.ax.xaxis.set_ticks_position('bottom')
cbar.ax.xaxis.set_label_position('bottom')
cbar.set_label('$U$  [m/s]', fontsize=14, labelpad=14)
cbar.update_ticks()
figname= "fig02.png"
plt.savefig(figpath+figname, dpi=150, bbox_inches='tight')
