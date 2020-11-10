import numpy as np
from netCDF4 import Dataset

import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib

import os

matplotlib.rcParams['xtick.labelsize'] = 12
matplotlib.rcParams['ytick.labelsize'] = 12
matplotlib.rcParams['axes.linewidth'] = 0.8
matplotlib.rcParams['axes.edgecolor']='black'
matplotlib.rcParams['font.family'] = 'serif'
# This requires a working LaTeX installation. If you don't have that,
# you may comment this line
matplotlib.rcParams['text.usetex']=True

figpath = "../figs/"

period = 10.3
slopes = [1.66, 2.0, 2.5, 3.0]
divs = [0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.0]
st = {'1.66':'5/3', '2.0':'2.0', '2.5':'2.5', '3.0':'3.0'}
bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=1, alpha=.8)

######################################################################
fig, axes = plt.subplots(1, 4, figsize=(18,4), sharey=True)
for i, slope in enumerate(slopes):
    ax = axes[i] 
    f0 = "K%sA0.0T10.3S0.01_hs.nc" %slope
    data0 = os.path.join('../data/synthetic/standard/gridded', f0)
    nc0 = Dataset(data0, 'r')
    hs0 = nc0.variables['hs'][:20]
    r = []
    for div in divs:
        if div not in [0.9, 0.95]:
            path = "../data/synthetic/standard/gridded/"
        else:
            path = "../data/synthetic/alpha9/gridded/"
        corr = []
        filename = "K%sA%sT10.3S0.01_hs.nc" %(slope, div)
        data = os.path.join(path, filename)
        nc = Dataset(data, 'r')
        hs = nc.variables['hs'][:20]
        corr = np.corrcoef(hs0.flatten(), hs.flatten())[0,1]
        r.append(corr)
    r = np.array(r)
    ax.plot(divs, r, '-o',lw=.5, color='steelblue')
    ax.grid(True, ls='dotted')
    ax.set_xticks(divs)
    ax.set_xticklabels(['0', '0.2', '0.4', '0.6', '0.8', '0.9', '', '1'])
    ax.text(0.1, 0.1, '$q^{-%s}$'%st['%s'%slope], fontsize=14, bbox=bbox_props)
fig.text(0.085, 0.45,'correlation coefficient', ha='center', va='center',
        rotation='vertical', fontsize=16)
fig.text(0.53, 0.0005,r'divergence fraction ($\alpha$)', ha='center', va='center', fontsize=16)
plt.subplots_adjust(wspace=0.03, hspace=0.0)
figname= "fig05.png" 
plt.savefig(figpath+figname, dpi=150, bbox_inches='tight')
