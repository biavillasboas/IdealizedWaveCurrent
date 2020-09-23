import numpy as np
from numpy import pi
from scipy.special import gammainc
import xarray as xr
import xrft
import dask

def calc_spec(phi, d1, d2):
    # wavenumber one (equals to dk1 and dk2)
    # time is the first dimension

    if phi.ndim == 2:
        n1, n2 = phi.shape
    else:
        nt, n1, n2 = phi.shape

    k = 2*pi*np.fft.rfftfreq(n2, d2)
    l = 2*pi*np.fft.fftshift(np.fft.fftfreq(n1, d1))
    kk, ll = np.meshgrid(k, l)
    dk = k[1] - k[0]
    dl = l[1] - l[0]

    kappa2 = kk**2 + ll**2
    kappa = np.sqrt(kappa2)

    phih = np.fft.rfft2(phi,axes=(-2,-1))
    spec = 2.*(phih*phih.conj()).real/ (dk*dl)\
                / (n1*n2)**2

    var_dens = np.fft.fftshift(spec.copy(),axes=-2)
    var_dens[:,:,0], var_dens[:,:,-1] = var_dens[:,:,0]/2., var_dens[:,:,-1]/2.

    return spec, var_dens, k, l, dk, dl


def calc_ispec(k, l, E):
    """ Calculates the azimuthally-averaged spectrum

        Parameters
        ===========
        - E is the two-dimensional spectrum
        - k is the wavenumber is the x-direction
        - l is the wavenumber in the y-direction

        Output
        ==========
        - kr: the radial wavenumber
        - Er: the azimuthally-averaged spectrum """
    
    dk = np.abs(k[2]-k[1])
    dl = np.abs(l[2]-l[1])

    k, l = np.meshgrid(k,l)

    wv = np.sqrt(k**2+l**2)

    if k.max()>l.max():
        kmax = l.max()
    else:
        kmax = k.max()

    dkr = np.sqrt(dk**2 + dl**2)
    kr =  np.arange(dkr/2.,kmax+dkr/2.,dkr)
    if E.ndim == 3:
        nt = E.shape[0]
        Er = np.zeros((nt, kr.size))
    else:
        Er = np.zeros((kr.size))

    for i in range(kr.size):
        fkr =  (wv>=kr[i]-dkr/2) & (wv<=kr[i]+dkr/2)
        dth = np.pi / (fkr.sum()-1)
        if E.ndim == 3:
            Er[:,i] = (E[:,fkr]*(wv[fkr])*dth).sum(axis=(-1))
        else:
            Er[i] = (E[fkr]*(wv[fkr])*dth).sum(axis=(-1))

    return kr, Er.squeeze()

def spectral_slope(k,E,kmin,kmax,stdE):
    ''' compute spectral slope in log space in
        a wavenumber subrange [kmin,kmax],
        m: spectral slope; mm: uncertainty'''

    fr = np.where((k>=kmin)&(k<=kmax))

    ki = np.matrix((np.log10(k[fr]))).T
    Ei = np.matrix(np.log10(np.real(E[fr]))).T
    dd = np.matrix(np.eye(ki.size)*((np.abs(np.log10(stdE)))**2))

    G = np.matrix(np.append(np.ones((ki.size,1)),ki,axis=1))
    Gg = ((G.T*G).I)*G.T
    m = Gg*Ei
    mm = np.sqrt(np.array(Gg*dd*Gg.T)[1,1])
    yfit = np.array(G*m)
    m = np.array(m)[1]

    return m, mm

def smith2007_decomposition(u, v, dx, dy):

    """ decompose the vector field (u,v) into potential (u_phi,v_phi)
        and rotational (u_psi,v_psi) fields using 2D FT following
         Smith (2008)
    !!!!!!!!!!!!!!! Doesn't work for odd kNyq, I have to fix that !!!!!!!!
    """

    ny,nx = u.shape
    k = np.fft.fftfreq(nx, dx)
    k = np.fft.fftshift(k)
    l = np.fft.fftfreq(ny, dy)
    l = np.fft.fftshift(l)

    kk,ll = np.meshgrid(k, l)
    THETA = (np.arctan2(ll, kk))

    U = np.fft.fft2(u,axes=(0, 1))
    V = np.fft.fft2(v,axes=(0, 1))

    P = U*np.cos(THETA) + V*np.sin(THETA)
    S = -U*np.sin(THETA) + V*np.cos(THETA)

    # back to physical space
    u_phi = np.real(np.fft.ifft2(P*np.cos(THETA),axes=(0,1)))
    v_phi = np.real(np.fft.ifft2(P*np.sin(THETA),axes=(0,1)))

    u_psi = np.real(np.fft.ifft2(-S*np.sin(THETA),axes=(0,1)))
    v_psi = np.real(np.fft.ifft2(S*np.cos(THETA),axes=(0,1)))

    assert np.allclose(u, (u_phi+u_psi), rtol=1e-5, atol=1e-5), "Decomposition error: rotational\
        and divergent u-components are not orthogonal"
    assert np.allclose(v, (v_phi+v_psi), rtol=1e-5, atol=1e-5), "Decomposition error: rotational\
        and divergent v-components are not orthogonal"

    return u_phi, v_phi, u_psi, v_psi


def helm_decompose_flow(u, v, dx):

    u_phi = u.copy()
    v_phi = u.copy() 
    u_psi = u.copy()
    v_psi = u.copy()

    if len(u.shape) == 3:
        nt = u.shape[0]
        for i in range(nt):
            u_phi[i], v_phi[i], u_psi[i], v_psi[i] = smith2007_decomposition(u[i],v[i], dx, dx)
    else:
        u_phi.values, v_phi.values, u_psi.values, v_psi.values = smith2007_decomposition(u, v, dx, dx)

    return u_phi, v_phi, u_psi, v_psi

def yNlu(sn, yN, ci):
    """ compute yN[l] yN[u], that is, the lower and
                upper limit of yN """

    # cdf of chi^2 dist. with 2*sn DOF
    cdf = gammainc(sn, sn*yN)

    # indices that delimit the wedge of the conf. interval
    fl = np.abs(cdf - ci).argmin()
    fu = np.abs(cdf - 1. + ci).argmin()

    return yN[fl], yN[fu]


def spec_error(E, sn, ci=0.95):
    dbin = 0.005
    yN = np.arange(0,2.+dbin,dbin)
    El, Eu = np.empty_like(E), np.empty_like(E)

    try:
        n = sn.size
    except AttributeError:
        n = 0
    if n:
        assert n == E.size, " *** sn has different size than E "

        for i in range(n):
            yNl, yNu = yNlu(sn[i], yN=yN, ci=ci)
            El[i] = E[i]/yNl
            Eu[i] = E[i]/yNu

    else:
        yNl, yNu = yNlu(sn, yN=yN, ci=ci)
        El = E/yNl
        Eu = E/yNu

    return El, Eu

def compute_spectra_ll_month(yymm):

    dsc = xr.open_dataset('../data/llc4320/gridded/llc_T10.3_20%s_cur.nc' %yymm)
    dsh = xr.open_dataset('../data/llc4320/gridded/llc_T10.3_20%s_hs.nc' %yymm)
    ni = np.where(dsc.x==200000)[0][0]
    nf = np.where(dsc.x==600000)[0][0]
    u_mix = dsc.ucur[:, ni:nf, ni:nf]
    v_mix = dsc.vcur[:, ni:nf, ni:nf]
    hs = dsh.hs[:, ni:nf, ni:nf]
    dx = (dsh.x.values[1] - dsh.x.values[0])
    u_phi, v_phi, u_psi, v_psi = helm_decompose_flow(u_mix, v_mix, dx)

    us_hat2 = xrft.power_spectrum(u_psi, dim=['x','y'], detrend='linear', window=True)
    vs_hat2 = xrft.power_spectrum(v_psi, dim=['x','y'], detrend='linear', window=True)
    ekes_hat2 = 0.5*(us_hat2 + vs_hat2)
    nk = us_hat2.shape[-1]
    ekes_y = ekes_hat2.sum(dim='freq_x')[:,nk//2:]*ekes_hat2.freq_x_spacing

    up_hat2 = xrft.power_spectrum(u_phi, dim=['x','y'], detrend='linear', window=True)
    vp_hat2 = xrft.power_spectrum(v_phi, dim=['x','y'], detrend='linear', window=True)
    ekep_hat2 = 0.5*(up_hat2 + vp_hat2)
    ekep_y = ekep_hat2.sum(dim='freq_x')[:,nk//2:]*ekep_hat2.freq_x_spacing

    hs_hat2 = xrft.power_spectrum(hs, dim=['x','y'], detrend='linear', window=True)
    nk = hs_hat2.shape[-1]
    hsy = hs_hat2.sum(dim='freq_x')[:,nk//2:]* hs_hat2.freq_x_spacing
    N = hsy.time.size
    ky = hsy.freq_y*1e3

    Es = ekes_y.mean(dim='time')*1e-3
    Ep = ekep_y.mean(dim='time')*1e-3
    Ehs = hsy.mean(dim='time')*1e-3
    Es, Ep, Ehs = dask.compute(Es, Ep, Ehs)
    return N, ky, Es, Ep, Ehs
