import numpy as np
import xarray as xr
from scipy.stats import linregress

import matplotlib.pyplot as plt

import pandas as pd
import seaborn as sns
import seawater as sw

import matplotlib
matplotlib.rcParams['xtick.labelsize'] = 22
matplotlib.rcParams['ytick.labelsize'] = 22
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['text.latex.preamble'] = [
    r'\usepackage{amsmath}',
    r'\usepackage{amssymb}']


ds1 = xr.open_dataset('../data/model_stats/S0.01_gridded_stats.nc')
ds1 = ds1.drop('index')
ds1.set_coords(['alpha', 'slope', 'seed', 'period'])
df1 = ds1.to_dataframe()

cg1 = 9.8*df1.t0m1_mean/4/np.pi
f = sw.f(31.5)
ind1=df1.alpha<0.9
a, b, __, __, __ = linregress(df1.vorticity[ind1]/f, cg1[ind1]*df1.hs_grad[ind1]/df1.hs_mean[ind1]/df1.slope[ind1]/f)


cmap = sns.light_palette("navy", 4, as_cmap=True)
markers = ['*', 'v', 's']

cs = plt.scatter(df1.vorticity/f, df1.hs_grad/df1.hs_mean/f, c = df1.slope, vmin=1.6, vmax=3, cmap=cmap)

fig, axes = plt.subplots(1, 3, figsize=(22,6))
ax1, ax2, ax3 = axes
for i, period in enumerate(df1.period.unique()):
    m = markers[i]
    d = df1[(ind1) & (df1.period==period)]
    cg = cg1[(ind1) & (df1.period==period)]
    ax1.scatter(d.vorticity/f, d.hs_grad/d.hs_mean/f, c = d.slope, alpha=.4, marker=m,
            vmin=1.6, vmax=3, label='$c=%.1f$ m/s'%(period*9.8/4/np.pi), cmap=cmap)
    ax2.scatter(d.vorticity/f, cg*d.hs_grad/d.hs_mean/f, c = d.slope, alpha=.4, marker=m,
            vmin=1.6, vmax=3, label='$c=%.1f$ m/s'%(period*9.8/4/np.pi), cmap=cmap)
    ax3.scatter(d.vorticity/f, cg*d.hs_grad/d.hs_mean/d.slope/f, c=d.slope, alpha=.4, marker=m,
            vmin=1.6, vmax=3, label='$c=%.1f$ m/s'%(period*9.8/4/np.pi), cmap=cmap)
ax3.plot(np.linspace(0, 0.3, 30), a*np.linspace(0, 0.3, 30) + b, '--', color='k')
ax1.set_ylim([1e-3, 0.11])
ax1.set_xlim([0, 0.35])
ax2.set_xlim([0, 0.35])
ax3.set_xlim([0, 0.35])
ax1.set_ylabel(r'$|\nabla H_s|_{rms} / \langle H_s \rangle f$', fontsize=26, labelpad=20)
ax2.set_ylabel(r'$c|\nabla H_s|_{rms} / \langle H_s\rangle f$', fontsize=26, labelpad=20)
ax3.set_ylabel(r'$c|\nabla H_s|_{rms} / \langle H_s\rangle \mathbb{S} f$', fontsize=26, labelpad=20)

leg = ax1.legend(loc=2, fontsize=22)
for h in leg.legendHandles:
    h.set_color('steelblue')
#    h.set_facecolor('w')
leg = ax2.legend(loc=2, fontsize=22)
for h in leg.legendHandles:
    h.set_color('steelblue')
    h.set_facecolor('w')
leg = ax3.legend(loc=2, fontsize=22)
for h in leg.legendHandles:
    h.set_color('steelblue')
    h.set_facecolor('w')
fig.text(0.51, 0.005, r'$\zeta_{rms}/f$', ha='center', va='center', fontsize=26)
# rect = [left, bottom, width, height]
cax1 = plt.axes([0.2, 0.25, 0.1, 0.015])
cbar1 = plt.colorbar(cs, cax=cax1, orientation='horizontal', ticks=[1.66, 2.0, 2.5, 3.0], format = '%.1f')
cbar1.ax.xaxis.set_ticks_position('bottom')
cbar1.ax.xaxis.set_label_position('bottom')
cbar1.set_label('$\mathbb{S}$',fontsize=22)
cbar1.ax.tick_params(labelsize=22)
cax2 = plt.axes([0.49, 0.25, 0.1, 0.015])
cbar2 = plt.colorbar(cs, cax=cax2, orientation='horizontal', ticks=[1.66, 2.0, 2.5, 3.0], format = '%.1f')
cbar2.ax.xaxis.set_ticks_position('bottom')
cbar2.ax.xaxis.set_label_position('bottom')
cbar2.set_label('$\mathbb{S}$',fontsize=22)
cbar2.ax.tick_params(labelsize=22)
cax3 = plt.axes([0.78, 0.25, 0.1, 0.015])
cbar3 = plt.colorbar(cs, cax=cax3, orientation='horizontal', ticks=[1.66, 2.0, 2.5, 3.0], format = '%.1f')
cbar3.ax.xaxis.set_ticks_position('bottom')
cbar3.ax.xaxis.set_label_position('bottom')
cbar3.set_label('$\mathbb{S}$',fontsize=22)
cbar3.ax.tick_params(labelsize=22)
cbar2.ax.set_xticklabels(['5/3', '2', '2.5', '3'])
cbar1.ax.set_xticklabels(['5/3', '2', '2.5', '3'])
cbar3.ax.set_xticklabels(['5/3', '2', '2.5', '3'])
fig.subplots_adjust(wspace=0.4, hspace=1.5)
plt.savefig('../figs/figS1.png',dpi=150, bbox_inches='tight')
