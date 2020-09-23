import os
import glob
import numpy as np
import pickle as pk
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm

matplotlib.rcParams['xtick.labelsize'] = 14
matplotlib.rcParams['ytick.labelsize'] = 14
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['text.usetex']=True

figpath = '../figs/'
datapath = '../data/model_stats/inverse'

nd = 10
months = ['%.2d' %d for d in range(1, 13)] 
corr = np.zeros([len(months), 250-nd])
for i, m in enumerate(months):
    finverse = open(glob.glob(os.path.join(datapath,
        'llc_inverse_model_llc_T10.3_20*%s.pkl'%m))[0], 'rb')
    ftrue = open(glob.glob(os.path.join(datapath,
        'llc_true_gradient_llc_T10.3_20*%s.pkl'%m))[0], 'rb')

    inverse_dict = pk.load(finverse)
    true_dict = pk.load(ftrue)
    inverse = inverse_dict[list(inverse_dict.keys())[0]][:, nd:-nd, :-nd]
    true = true_dict[list(true_dict.keys())[0]][:, nd:-nd, :-nd]

    _, _, Nx = true.shape
    for nx in range(Nx):
        corr[i, nx] = np.corrcoef(inverse[:, :, nx].flatten(),
                true[:, :, nx].flatten())[0, 1]
    print('%s'%m)
    

t = np.arange(1,13)
tlabels = ['Jan', 'Mar', 'May', 'Jul', 'Sep', 'Nov']
x = np.arange(Nx)*2.5

fig, ax = plt.subplots(1, 1, figsize=(8,6))
cs = ax.pcolormesh(x, months, corr, vmin=0.35, vmax=0.75, cmap=cm.magma)
ax.set_xlabel('x distance [km]', fontsize=16, labelpad=14)
ax.set_yticks([1, 3, 5, 7, 9, 11])
ax.set_yticklabels(tlabels, rotation=45)
# rect = [left, bottom, width, height]
cax = plt.axes([0.92, 0.2, 0.02, 0.6])
cbar = plt.colorbar(cs, cax=cax, extend='both', format = '%.2f', ticks=[0.35, 0.45, 0.55, 0.65, 0.75, ])
cbar.set_label('correlation coefficient', fontsize=14,  labelpad=14)
cbar.ax.tick_params(labelsize=14)
plt.savefig(figpath+'figS2.png', dpi=150, bbox_inches='tight')
