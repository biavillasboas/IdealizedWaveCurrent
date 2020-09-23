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
######################################################################
hs_lin = np.linspace(0.8, 1.2, 30)
bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=1, alpha=.8)
std = 0.01
st = {'1.66':'5/3', '2.0':'2.0', '2.5':'2.5', '3.0':'3.0'}
seed_idx =13
i=1
fig = plt.figure(figsize=(9.9,10))
for slope in slopes:
    for div in divs:

        filename = "K%sA%sT%sS%s_hs.nc" %(slope, div, period, std)
        data = os.path.join(path, filename)
        nc = Dataset(data, 'r')

        hs = nc.variables['hs'][seed_idx]
        x = (nc.variables['x'][:] - nc.variables['x'][0])*1e-3
        y = (nc.variables['y'][:] - nc.variables['y'][0])*1e-3

        ax = plt.subplot(4,4,i) 
        cs = ax.contourf(x, y, hs, hs_lin, extend='both', cmap='cmo.curl')

        ax.set_xticks([50,250, 450])
        ax.set_yticks([50,250, 450])
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.tick_params(labelbottom='off',labeltop='on')

        if i not in [1, 2, 3, 4]:
            ax.set_xticklabels([])
            ax.set_xticks([])
            plt.setp(ax.get_xticklabels(), visible=False)

        if i not in [1,5,9,13]:
            ax.set_yticklabels([])
            ax.set_yticks([])
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
cbar = plt.colorbar(cs,cax=cax, extend='both', orientation='horizontal',format = '%.1f')
cbar.ax.xaxis.set_ticks_position('bottom')
cbar.ax.xaxis.set_label_position('bottom')
cbar.set_label('$H_s$  [m]', fontsize=14, labelpad=14)
tick_locator = ticker.MaxNLocator(nbins=4)
cbar.locator = tick_locator
cbar.update_ticks()
figname= "fig04.png"
plt.savefig(figpath+figname, dpi=150, bbox_inches='tight')
