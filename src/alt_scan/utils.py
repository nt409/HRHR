import numpy as np
from math import floor

from model.params import PARAMS
from model.utils import get_rfd



def get_alt_scan_params(n_its, index):
    rfs = np.logspace(-10, -4, n_its)
    omegas = np.linspace(0.4, 1, n_its)
    thetas = np.linspace(4,  12, n_its)


    ii = int( floor(index/(n_its**2)) )
    jj = int( floor((index % 25)/(n_its)) )
    kk = int( index % n_its )

    rf1 = rfs[ii]
    rf2 = 1e-5
    om1 = omegas[jj]
    om2 = 1
    thet1 = thetas[kk]
    thet2 = 9.6

    rfd = get_rfd(rf1, rf2)
    
    primary_inoc = dict(RR=rfd, RS=rf1, SR=rf2, SS=1-rf1-rf2-rfd)

    fcide_parms = dict(
        omega_1 = om1,
        omega_2 = om2,
        theta_1 = thet1,
        theta_2 = thet2,
        delta_1 = PARAMS.delta_1,
        delta_2 = PARAMS.delta_2,
        )
    
    return primary_inoc, fcide_parms, ii, jj, kk