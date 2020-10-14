import numpy as np
import xarray as xr
from scipy.stats import linregress

import matplotlib.pyplot as plt

import pandas as pd
import seaborn as sns
import seawater as sw

import matplotlib
matplotlib.rcParams['xtick.labelsize'] = 14
matplotlib.rcParams['ytick.labelsize'] = 14
matplotlib.rcParams['font.family'] = 'serif'
# This requires a working LaTeX installation. If you don't have that,
# you may comment the lines below and change the plot ylabel
matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['text.latex.preamble'] = [
            r'\usepackage{amsmath}',
                r'\usepackage{amssymb}']

figpath = '../figs/'

dsl = xr.open_dataset('../data/model_stats/llc_gridded_stats.nc')
dfl = dsl.to_dataframe()

ds1 = xr.open_dataset('../data/model_stats/S0.01_gridded_stats.nc')
ds1 = ds1.drop('index')
ds1.set_coords(['alpha', 'slope', 'seed', 'period'])
df1 = ds1.to_dataframe()

ds5 = xr.open_dataset('../data/model_stats/S0.005_gridded_stats.nc')
ds5 = ds5.drop('index')
ds5.set_coords(['alpha', 'slope', 'seed', 'period'])
df5 = ds5.to_dataframe()
df = pd.concat([df1, df5])

cg = 9.8*df.t0m1_mean/4/np.pi
cg1 = 9.8*df1.t0m1_mean/4/np.pi
cg5 = 9.8*df5.t0m1_mean/4/np.pi
c_llc = 9.8*dfl.t0m1_mean/4/np.pi
f = sw.f(31.5)

ind=df.alpha<0.9
ind1=df1.alpha<0.9
ind5=df5.alpha<0.9

corr_df = np.corrcoef(df.vorticity[ind].values, cg[ind]*df.hs_grad[ind]/df.hs_mean[ind]/df.slope[ind])[0, 1]
cor_llc = np.corrcoef(dfl.vorticity, c_llc*dfl.hs_grad/dfl.hs_mean/dfl.EKE_psi_slope)[0,1]
a, b, __, __, __ = linregress(df.vorticity[ind]/f, cg[ind]*df.hs_grad[ind]/df.hs_mean[ind]/df.slope[ind]/f)

############################################
# Normalized scatterplot
############################################
plt.figure(figsize=(8,8))
plt.plot(df1.vorticity[ind1]/f, cg1[ind1]*df1.hs_grad[ind1]/df1.hs_mean[ind1]/df1.slope[ind1]/f, 'o',
                        color='steelblue', alpha=.5, label='KE = 0.01 m$^2$/s$^2$')
plt.plot(df5.vorticity[ind5]/f, cg5[ind5]*df5.hs_grad[ind5]/df5.hs_mean[ind5]/df5.slope[ind5]/f, 'o',
                        color='darkblue', alpha=.4, label='KE = 0.005 m$^2$/s$^2$')
plt.plot(dfl.vorticity/f, c_llc*dfl.hs_grad/dfl.hs_mean/dfl.EKE_psi_slope/f, '+',
                        alpha=0.4, color='k', label='LLC4320')
plt.plot(np.linspace(0, 0.3, 30), a*np.linspace(0, 0.3, 30) + b, '--', color='k')
plt.legend(facecolor='w', fontsize=14)
plt.grid(ls='dotted')
plt.xlabel('$\zeta_{rms}/f$', fontsize=18)
plt.xlim([0, 0.3])
plt.ylim([0, 0.4])
plt.ylabel(r'$c|\nabla H_s|_{rms} / \langle H_s\rangle \mathbb{S} f$', fontsize=18, labelpad=20)
plt.savefig(figpath+'fig12.png', dpi=150, bbox_inches='tight')
