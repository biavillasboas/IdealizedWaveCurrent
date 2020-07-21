import os
import numpy as np
import scipy as sp

import datetime
import netCDF4
from netCDF4 import Dataset, date2num
from numpy import pi



def create_synthetic_flow_nc(seed, slope, alpha, L, dx, kmin, kmax, mean_eke, output_path):
    """ Creates synthetic flow and saves as a NetCDF file.

    This function produces three netCDF files containing the u and v 
    components of the synthetic velocity. The files are saved with the
    prefixes u_mix, u_psi, and u_phi, corresponding to:

    u_mix: Flow with prescribed divergence fraction `alpha` (see below)
    u_psi: Only the rotational component of u_mix
    u_phi: Only the potential component of u_mix

    file names follow the pattern: 'u_%s_k%s_alpha%s_seed%d.nc'
    %(prefix, slope, alpha, seed)

    Parameters
    ----------
    seed : int
	Number to be used as a seed for numpy.random.RandomState(seed) to generate
	the phase of the synthetic flow.
    slope : float
	Slope of the isotropic kinetic energy (KE) spectrum of the synthetic flow.
	The resulting spectrum will follow a $k^-slope$ power law, so for a red
	KE spectrum (as usually observed in the ocean) `slope` should be positive.
    alpha : float
	Fraction of total kinetic energy in the divergent (potential) component of
	the flow (between 0 and 1). If `alpha` = 0, the flow has no divergence 
	(it is purely rotational) if `alpha`= 1 the flow has no vorticity (purely
	potential).
    L : float
	Length of the square domain in meters.
    dx : float
	Regular grid spacing in meters.
    kmin : float
	Lowest wavenumber containing energy (in radians per meter). Note that this
	will only affect the result if `kmin` is higher than 2pi/L.
    kmax : float
	Highest wavenumber containing energy (in radians per meter). Note that this
	will only affect the result if `kmax` is lower than the Nyquest wavenuber.
    mean_eke : float
	Mean kinetic energy of the resulting synthetic flow (in m^s/s^2).
    output_path : string
	Path where to save the output netCDF file.
    """

    output_name_psi = 'u_%s_k%s_alpha%s_seed%d.nc' %('psi', slope, alpha, seed)
    output_psi = os.path.join(output_path, output_name_psi)
    output_name_phi = 'u_%s_k%s_alpha%s_seed%d.nc' %('phi', slope, alpha, seed)
    output_phi = os.path.join(output_path, output_name_phi)
    output_name_mix = 'u_%s_k%s_alpha%s_seed%d.nc' %('mix', slope, alpha, seed)
    output_mix = os.path.join(output_path, output_name_mix)

    x = np.arange(0, L, dx)
    y = x

    us, vs, up, vp = potential_solenoidal(L, dx, slope, kmin, kmax, seed)
    u_psi, v_psi, u_phi, v_phi, u_mix, v_mix = synthetic_flow(us, vs,
                                                              up, vp, mean_eke, alpha)
    save_netcdf_fields(x=x, y=y, u=u_mix, v=v_mix, time=None, seed=seed, slope=slope,
                       alpha=alpha, flow_type='mix', output=output_mix)
    save_netcdf_fields(x=x, y=y, u=u_psi, v=v_psi, time=None, seed=seed, slope=slope,
                       alpha=alpha, flow_type='psi', output=output_psi)
    save_netcdf_fields(x=x, y=y, u=u_phi, v=v_phi, time=None, seed=seed, slope=slope,
                       alpha=alpha, flow_type='phi', output=output_phi)

    return


def potential_solenoidal(L, dx, slope, kmin, kmax, seed=10):
    """ Creates purely rotational and purely divergent velocity fields.
    
    Returns four arrays corresponding to the u and v components of a purely 
    rotational and a purely divergent velocity field with variance 
    normalized to 1 and prescribed kinetic energy spectral slope `slope`. The
    rotational and divergent flows are computed from the same function, used 
    as a stream function (for the rotational flow) and a velocity potential
    for the divergent component.

    Parameters
    ----------
    L : float
	Length of the square domain in meters.
    dx : float
        Regular grid spacing in meters.
    slope : float
        Slope of the isotropic kinetic energy (KE) spectrum of the synthetic flow.
	The resulting spectrum will follow a $k^-slope$ power law, so for a red
	KE spectrum (as usually observed in the ocean) `slope` should be positive.
    kmin : float
	Lowest wavenumber containing energy (in radians per meter). Note that this
	will only affect the result if `kmin` is higher than 2pi/L.
    kmax : float
	Highest wavenumber containing energy (in radians per meter). Note that this
	will only affect the result if `kmax` is lower than the Nyquest wavenuber.
    seed : int
	Number to be used as a seed for numpy.random.RandomState(seed) to generate
	the phase of the synthetic flow.

    Returns
    -------
    u_psi : ndaray
        u-component of the purely rotational velocity
    v_psi : ndaray
        v-component of the purely rotational velocity
    u_phi : ndaray
        u-component of the purely divergent velocity
    v_phi : ndaray
        v-component of the purely divergent velocity

    Notes
    -----
    The kinetic energy 0.5*(u^2 + v^2) of both the rotational (psi) and
    divergent (phi) components is normalized to one i.e., 
    0.5*(u_psi^2 + v_psi^2) = 0.5*(u_phi^2 + v_phi^2) = 1

    """

    N =  int(L/dx) # number of points
    l = 2*pi*np.fft.fftshift(np.fft.fftfreq(N, dx))
    k = l
    kk, ll = np.meshgrid(k, l)
    K = (kk**2 + ll**2)**.5 # isotropic wavenumber
    klim = k.max()
    K = np.ma.masked_array(K, K<kmin)
    K = np.ma.masked_array(K, K>kmax)
    S0 =(1/K**slope)/K # extra factor of K from polar Jacobian
    Nk, Nl = len(k), len(l)

    rand_seed = np.random.RandomState(seed)
    pha = 2*np.pi*(2*rand_seed.rand(Nl, Nk) - 1)
    F_hat = (np.cos(pha) + 1j*np.sin(pha))*S0**.5
    F_hat[K>klim] = 0. + 0.j

    Fy_hat = np.fft.fftshift(F_hat * 1.0j*ll/K)
    Fx_hat = np.fft.fftshift(F_hat* 1.0j*kk/K)
    u_psi = np.fft.ifft2(-Fy_hat).real
    v_psi = np.fft.ifft2(Fx_hat).real
    u_phi = np.fft.ifft2(Fx_hat).real
    v_phi = np.fft.ifft2(Fy_hat).real

    phi_var = (u_phi.var() + v_phi.var())**.5
    u_phi = u_phi/phi_var
    v_phi = v_phi/phi_var

    psi_var = (u_psi.var() + v_psi.var())**.5
    u_psi = u_psi/psi_var
    v_psi = v_psi/psi_var

    return u_psi, v_psi, u_phi, v_phi


def synthetic_flow(u_psi, v_psi, u_phi, v_phi, mean_eke, alpha):
    """Combines rotational and divergent velocities

    Takes u and v components of purely rotational and purely divergent
    velocities and combine them using a prescribed divergence fraction
    `alpha` and normalizes the total kinetic energy to `mean_eke`.
    Returns the u and v components of the rotational, divergent and 
    combined velocity.

    Parameters
    ----------
    u_psi : ndarray
        u-component of the rotational velocity.
    v_psi : ndarray
        v-component of the rotational velocity.
    u_phi : ndarray
        u-component of the divergent velocity.
    v_phi : ndarray
        v-component of the divergent velocity.
    mean_eke : float
        Mean kinetic energy of the resulting synthetic flow (in m^s/s^2).
    alpha : float
        Fraction of total kinetic energy in the divergent (potential) component of
        the flow (between 0 and 1). If `alpha` = 0, the flow has no divergence
        (it is purely rotational) if `alpha`= 1 the flow has no vorticity (purely
        potential).

    Returns
    -------
    u_psi : ndaray
        u-component of the purely rotational velocity
    v_psi : ndaray
        v-component of the purely rotational velocity
    u_phi : ndaray
        u-component of the purely divergent velocity
    v_phi : ndaray
        v-component of the purely divergent velocity
    u_mix : ndaray
        u-component of the total velocity
    v_mix : ndaray
        v-component of the total velocity

    Notes
    -----

    The input u_psi, v_psi, u_phi, v_phi have their variance normalized to 1,
    such that  0.5*(u_psi^2 + v_psi^2) = 0.5*(u_phi^2 + v_phi^2) = 1, while
    the output are normalized such that
    
    mean_eke = 
    `alpha`*0.5*(u_phi^2 + v_phi^2) + (1-`alpha`)*0.5*(u_psi^2 + v_psi^2)

    """

    norm = (2*mean_eke)**.5

    u_psi = ((1-alpha)**.5)*norm*u_psi
    v_psi = ((1-alpha)**.5)*norm*v_psi

    u_phi = (alpha**.5)*norm*u_phi
    v_phi = (alpha**.5)*norm*v_phi

    EKE_phi = 0.5*(u_phi**2 + v_phi**2).mean()
    EKE_psi = 0.5*(u_psi**2 + v_psi**2).mean()

    assert np.allclose(EKE_phi, alpha*mean_eke), "Normalization error"
    assert np.allclose(EKE_psi, (1-alpha)*mean_eke), "Normalization error"

    u_mix = u_psi + u_phi
    v_mix = v_psi + v_phi
    EKE_mix = 0.5*(u_mix**2 + v_mix**2).mean()
    assert np.allclose(EKE_mix, mean_eke), "Normalization error"

    return u_psi, v_psi, u_phi, v_phi, u_mix, v_mix

