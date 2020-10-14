import sys
sys.path.append('../tools')

import glob
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib

import seawater as sw
import cmocean.cm as cmo

from inverse_tools import *

matplotlib.rcParams['axes.linewidth'] = 0.8
matplotlib.rcParams['axes.edgecolor']='black'
matplotlib.rcParams['font.family'] = 'serif'
# This requires a working LaTeX installation. If you don't have that,
# you may comment this line
matplotlib.rcParams['text.usetex']=True

figpath = '../figs/'

dir_file = '../data/llc4320/gridded/llc_T10.3_201110_dir.nc'
tm_file = '../data/llc4320/gridded/llc_T10.3_201110_t0m1.nc'
cur_file = '../data/llc4320/gridded/llc_T10.3_201110_cur.nc'

ncdir = Dataset(dir_file, 'r')
ww3_direction = ncdir.variables['dir'][0]
theta = np.radians(dir360_dir180(direction_from_to(az2trig(ww3_direction))))

nctm = Dataset(tm_file, 'r')
tm = nctm.variables['t0m1'][0]

nccur = Dataset(cur_file, 'r')
u = nccur.variables['ucur'][0]
x = nccur.variables['x'][:]*1e-3
y = nccur.variables['y'][:]*1e-3
v = nccur.variables['vcur'][0]
U = (u**2 + v**2)**.5

alpha = np.arctan2(v, u)
gamma = abs(alpha-theta)
Uk = U*np.cos(gamma)
Uky, Ukx = np.gradient(Uk)
thetay, thetax = np.gradient(theta)
gradk_theta = grad_k(theta, thetax, thetay)
gradn_Uk = grad_n(theta, Ukx, Uky)

k = (2*pi/tm)**2/9.8
cg = 9.8*tm/(4*np.pi)

model_theta = -cg*gradk_theta
##########################################
f = sw.f(31.5)
dx = 2.5
norm = 2.5e3*1e-4 # scale gradient by dx and  1e-4
bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=1, alpha=.9)

lin = np.linspace(-0.4, 0.4, 30)
fig, axes = plt.subplots(1,3,sharey='all', figsize=(18,6))
ax1 = axes.flatten()[0]
csu = ax1.contourf(x[:-10], y[:-10], U[:-10, :-10], np.linspace(0, 0.5), cmap='GnBu_r', extend='both')
cax = plt.axes([0.16, 0.05, 0.19, 0.02])
cbar = plt.colorbar(csu,cax=cax, extend='both', orientation='horizontal',format = '%.1f')
tick_locator = ticker.MaxNLocator(nbins=5)
cbar.locator = tick_locator
cbar.ax.tick_params(labelsize=18)
cbar.set_label('U  [m/s]', fontsize=16)
cbar.update_ticks()
ax1.xaxis.tick_top()
ax1.set_xticks([0,100, 200, 300, 400, 500])
ax1.set_xticklabels([0,100, 200, 300, 400, 500], fontsize=18)
ax1.set_yticks([0,100, 200, 300, 400, 500])
ax1.set_yticklabels([0,100, 200, 300, 400, 500], fontsize=18)
ax1.text(30, 540, 'A', fontsize=22, bbox=bbox_props)
ax2 = axes.flatten()[1]
ax2.contourf(x[:-10], y[:-10], gradn_Uk[:-10, :-10]/norm, lin, cmap=cmo.balance, extend='both')
ax2.xaxis.tick_top()
ax2.set_xticks([0,100, 200, 300, 400, 500])
ax2.set_xticklabels([0,100, 200, 300, 400, 500], fontsize=18)
ax2.text(465, 50, r'$\nabla_{\bot} U^k$', fontsize=22, bbox=bbox_props)
ax2.text(30, 540, 'B', fontsize=22, bbox=bbox_props)
ax3 = axes.flatten()[2]
cs = ax3.contourf(x[:-10], y[:-10], model_theta[:-10, :-10]/norm, lin, cmap=cmo.balance, extend='both')
cax2 = plt.axes([0.5, 0.05, 0.3, 0.02])
cbar2 = plt.colorbar(cs,cax=cax2, extend='both', orientation='horizontal',format = '%.1f')
tick_locator = ticker.MaxNLocator(nbins=5)
cbar2.locator = tick_locator
cbar2.ax.tick_params(labelsize=18)
cbar2.update_ticks()
cbar2.set_label(r'velocity gradient $\times 10^{4}$ [1/s]', fontsize=16)
ax3.text(30, 540, 'C', fontsize=22, bbox=bbox_props)
ax3.text(425, 50, r'$-c\nabla_{\parallel} \theta_w$', fontsize=22, bbox=bbox_props)
ax3.xaxis.tick_top()
ax3.set_xticks([0,100, 200, 300, 400, 500])
ax3.set_xticklabels([0,100, 200, 300, 400, 500], fontsize=18)
plt.subplots_adjust(wspace=0.05, hspace=0.0)
fig.text(0.065, 0.5, 'distance [km]', ha='center', va='center', rotation='vertical', fontsize=20)
fig.text(0.51, 0.98, 'distance [km]', ha='center', va='center', fontsize=20)
plt.savefig(figpath+'fig13.png', dpi=150, bbox_inches='tight')
