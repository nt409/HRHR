import numpy as np


from utils.functions import changing_fcide_dose_curve, changing_fcide_curve_asymp, \
    changing_fcide_sexp_asymp_curv

from utils.plotting.figures import fcide_grid

curv_dose = True
curv_asymp = True
sex_prop_asymp = True

NY = 30
N = 5

# plot

if curv_dose:
    # should use same N to emphasise symmetry
    doses = np.linspace(0, 1, N+1)[1:]
    curvatures = np.linspace(0, 10, N+1)[1:]
    rf = 10**(-5)

    x, y, z, filename_img = changing_fcide_dose_curve(doses, curvatures, rf, NY)
    labels = dict(x="Dose", y="Curvature", cbar="Failure year")
    fcide_grid(x, y, z, filename_img, labels)


if curv_asymp:
    # can use different N
    asymps = np.linspace(0, 1, N+1)[1:]
    curvatures = np.linspace(0, 10, N+1)[1:]
    rf = 10**(-5)

    x, y, z, filename_img = changing_fcide_curve_asymp(curvatures, asymps, rf, NY)
    labels = dict(x="Curvature", y="Asymptote", cbar="Failure year")
    fcide_grid(x, y, z, filename_img, labels)


if sex_prop_asymp:
    # can use different N
    asymps = np.linspace(0, 1, N+1)[1:]
    curvatures = np.linspace(0, 10, N+1)[1:]
    sex_props = np.linspace(0, 1, N)
    rf = 10**(-5)

    x, y, z, conf_str = changing_fcide_sexp_asymp_curv(sex_props,
                                            asymps, curvatures, rf, NY)
    labels = dict(x="Sexual proportion", y="Asymptote", cbar="Curvature")
    fcide_grid(x, y, z, conf_str, labels)

