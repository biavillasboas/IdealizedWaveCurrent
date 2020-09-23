import glob
import os
import sys
utils = os.path.expanduser('~/work/projects/python/utils')

import pickle
import numpy as np
from netCDF4 import Dataset
from inverse_tools import *

dir_files = sorted(glob.glob('../data/llc4320/gridded/llc_*dir.nc'))
tm_files = sorted(glob.glob('../data/llc4320/gridded/llc_*t0m1.nc'))
cur_files = sorted(glob.glob('../data/llc4320/gridded/llc_*cur.nc'))
N = len(dir_files)
dx = 2.5
figpath = "/home/avillasb/work/projects/research/synthetic_wavecurrent/figs/paper/llc/"
norm = 2.5e3*1e-4
bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=1, alpha=.9)

for n in range(N):
    inverse_model = {}
    true_gradient = {}
    ncdir = Dataset(dir_files[n], 'r')
    x = ncdir.variables['x'][:]*1e-3
    y = ncdir.variables['x'][:]*1e-3
    directions = ncdir.variables['dir'][:]
    nccur = Dataset(cur_files[n], 'r')
    us = nccur.variables['ucur'][:]
    vs = nccur.variables['vcur'][:]
    ns, ny, nx = directions.shape
    nctm = Dataset(tm_files[n], 'r')

    fname = os.path.splitext(os.path.basename(cur_files[n]))[0][:-4]
    inverse_model[fname] = np.zeros([ns, len(y), len(x)])
    true_gradient[fname] = np.zeros([ns, len(y), len(x)])
    
    for i in range(ns):

        u = us[i]
        v = us[i]
        U = (u**2 + v**2)**.5

        ww3_direction = directions[i]
        theta = np.radians(dir360_dir180(direction_from_to(az2trig(ww3_direction))))

        alpha = np.arctan2(v, u)
        gamma = abs(alpha-theta)
        Uk = U*np.cos(gamma)
        Uky, Ukx = np.gradient(Uk)
        thetay, thetax = np.gradient(theta)
        gradk_theta = grad_k(theta, thetax, thetay)
        gradn_Uk = grad_n(theta, Ukx, Uky)

        tm = nctm.variables['t0m1'][i]
        cg = 9.8*tm/(4*np.pi)

        inverse_theta = -cg*gradk_theta

        inverse_model[fname][i] =  inverse_theta
        true_gradient[fname][i] =  gradn_Uk


    output_inv = open("../data/model_stats/inverse/llc_inverse_model_%s.pkl" %fname, "wb")
    output_true = open("../data/model_stats/inverse/llc_true_gradient_%s.pkl" %fname, "wb")
    pickle.dump(inverse_model, output_inv)
    pickle.dump(true_gradient, output_true)
    output_inv.close()
    output_true.close()
    print('processed file %s' %fname)
    ncdir.close()
    nccur.close()
    nctm.close()
