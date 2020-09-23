import sys
sys.path.append('../tools')

import numpy as np

import xrft
import xarray as xr

import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['axes.linewidth'] = 0.8
matplotlib.rcParams['axes.edgecolor']='black'
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['text.usetex']=True

from spectra_tools import compute_spectra_ll_month, spec_error

figpath = '../figs/'


######################################
N, ky, Es_jan, Ep_jan, Ehs_jan = compute_spectra_ll_month('1201')
N, ky, Es_jul, Ep_jul, Ehs_jul = compute_spectra_ll_month('1207')

N, ky, Es_oct, Ep_oct, Ehs_oct = compute_spectra_ll_month('1110')
N, ky, Es_mar, Ep_mar, Ehs_mar = compute_spectra_ll_month('1203')

Es_jan_l, Es_jan_u = spec_error(Es_jan, N, ci=0.95)
Es_jul_l, Es_jul_u = spec_error(Es_jul, N, ci=0.95)
Es_oct_l, Es_oct_u = spec_error(Es_oct, N, ci=0.95)
Es_mar_l, Es_mar_u = spec_error(Es_mar, N, ci=0.95)

Ep_jan_l, Ep_jan_u = spec_error(Ep_jan, N, ci=0.95)
Ep_jul_l, Ep_jul_u = spec_error(Ep_jul, N, ci=0.95)
Ep_oct_l, Ep_oct_u = spec_error(Ep_oct, N, ci=0.95)
Ep_mar_l, Ep_mar_u = spec_error(Ep_mar, N, ci=0.95)

Ehs_jan_l, Ehs_jan_u = spec_error(Ehs_jan, N, ci=0.95)
Ehs_jul_l, Ehs_jul_u = spec_error(Ehs_jul, N, ci=0.95)
Ehs_oct_l, Ehs_oct_u = spec_error(Ehs_oct, N, ci=0.95)
Ehs_mar_l, Ehs_mar_u = spec_error(Ehs_mar, N, ci=0.95)

######################################
k1 = 1/200
k2 = 1/15
fig, axes = plt.subplots(2, 2, sharex='all', sharey='all', figsize=(8,12))
ax1, ax2, ax3, ax4 = axes.flatten()

ax1.loglog(ky, Es_jan, color='steelblue', ls='-', label=r'$\tilde E^\psi$ Jan')
ax1.fill_between(ky, Es_jan_l, Es_jan_u, color='steelblue', alpha=0.25)
ax1.loglog(ky, Es_jul, color='tomato', ls='-', label=r'$\tilde E^\psi$ Jul')
ax1.fill_between(ky, Es_jul_l, Es_jul_u, color='tomato', alpha=0.25)
hsjan, = ax1.loglog(ky, Ehs_jan/10, color='steelblue', ls='--')
ax1.fill_between(ky, Ehs_jan_l/10, Ehs_jan_u/10, color='steelblue', alpha=0.25)
hsjul, = ax1.loglog(ky, Ehs_jul/10, color='tomato', ls='--')
ax1.fill_between(ky, Ehs_jul_l/10, Ehs_jul_u/10, color='tomato', alpha=0.25)
ax1.axvspan(k1, k2, alpha=0.1, color='grey')
leg1 = ax1.legend(loc='upper right', framealpha=1, fontsize=13)
leg2 = ax1.legend([hsjan, hsjul],['$H_s$ Jan','$H_s$ Jul'], loc='lower left', framealpha=1, fontsize=13)
ax1.add_artist(leg1)

ax2.loglog(ky, Ep_jan, color='steelblue', ls='-', label=r'$\tilde E^\phi$ Jan')
ax2.fill_between(ky, Ep_jan_l, Ep_jan_u, color='steelblue', alpha=0.25)
ax2.loglog(ky, Ep_jul, color='tomato', ls='-', label=r'$\tilde E^\phi$ Jul')
ax2.fill_between(ky, Ep_jul_l, Ep_jul_u, color='tomato', alpha=0.25)
hsjan, = ax2.loglog(ky, Ehs_jan/10, color='steelblue', ls='--')
ax2.fill_between(ky, Ehs_jan_l/10, Ehs_jan_u/10, color='steelblue', alpha=0.25)
hsjul, = ax2.loglog(ky, Ehs_jul/10, color='tomato', ls='--')
ax2.fill_between(ky, Ehs_jul_l/10, Ehs_jul_u/10, color='tomato', alpha=0.25)
ax2.axvspan(k1, k2, alpha=0.1, color='grey')
leg1 = ax2.legend(loc='upper right', fontsize=13, framealpha=1)
leg2 = ax2.legend([hsjan, hsjul],['$H_s$ Jan','$H_s$ Jul'], loc='lower left', fontsize=13, framealpha=1)
ax2.add_artist(leg1)

ax3.loglog(ky, Es_oct, color='mediumslateblue', ls='-', label=r'$\tilde E^\psi$ Oct')
ax3.fill_between(ky, Es_oct_l, Es_oct_u, color='mediumslateblue', alpha=0.25)
ax3.loglog(ky, Es_mar, color='darkorange', ls='-', label=r'$\tilde E^\psi$ Mar')
ax3.fill_between(ky, Es_mar_l, Es_mar_u, color='darkorange', alpha=0.25)
hsoct, = ax3.loglog(ky, Ehs_oct/10, color='mediumslateblue', ls='--')
ax3.fill_between(ky, Ehs_oct_l/10, Ehs_oct_u/10, color='mediumslateblue', alpha=0.25)
hsmar, = ax3.loglog(ky, Ehs_mar/10, color='darkorange', ls='--')
ax3.fill_between(ky, Ehs_mar_l/10, Ehs_mar_u/10, color='darkorange', alpha=0.25)
ax3.axvspan(k1, k2, alpha=0.1, color='grey')
leg1 = ax3.legend(loc='upper right', fontsize=13, framealpha=1)
leg2 = ax3.legend([hsoct, hsmar],['$H_s$ Oct','$H_s$ Mar'], loc='lower left', fontsize=13, framealpha=1)
ax3.add_artist(leg1)

ax4.loglog(ky, Ep_oct, color='mediumslateblue', ls='-', label=r'$\tilde E^\phi$ Oct')
ax4.fill_between(ky, Ep_oct_l, Ep_oct_u, color='mediumslateblue', alpha=0.25)
ax4.loglog(ky, Ep_mar, color='darkorange', ls='-', label=r'$\tilde E^\phi$ Mar')
ax4.fill_between(ky, Ep_mar_l, Ep_mar_u, color='darkorange', alpha=0.25)
hsoct, = ax4.loglog(ky, Ehs_oct/10, color='mediumslateblue', ls='--')
ax4.fill_between(ky, Ehs_oct_l/10, Ehs_oct_u/10, color='mediumslateblue', alpha=0.25)
hsmar, = ax4.loglog(ky, Ehs_mar/10, color='darkorange', ls='--')
ax4.fill_between(ky, Ehs_mar_l/10, Ehs_mar_u/10, color='darkorange', alpha=0.25)
ax4.axvspan(k1, k2, alpha=0.1, color='grey')
leg1 = ax4.legend(loc='upper right', fontsize=13, framealpha=1)
leg2 = ax4.legend([hsoct, hsmar],['$H_s$ Oct','$H_s$ Mar'], loc='lower left', fontsize=13, framealpha=1)
ax4.add_artist(leg1)
for ax in axes.flatten():
    ax.grid(which='both', ls='dotted')
    ax.set_xlim([1e-3, 1e-1])
    ax.set_ylim([1e-6, 1e-1])
    ax.set_yticks([1e-6, 1e-4, 1e-2])
    ax.set_xticks([1e-2, 1e-1])
plt.subplots_adjust(wspace=0.1, hspace=0.1)
plt.annotate('Wavenumber [cycles/km]', xy=(0.38,0.03), xycoords="figure fraction", fontsize=12)
plt.annotate('PSD [m$^2$/cycle/km] or [m$^2$/s$^2$/cycle/km]', xy=(0.02,0.35), xycoords="figure fraction", fontsize=12, rotation =90)
plt.savefig(figpath+'fig11.png', bbox_inches='tight', dpi=150)
