import numpy as np

from utils.functions import changing_fcide_dose_curve, changing_fcide_curve_asymp
from utils.plotting import fcide_grid

curv_dose = False
curv_asymp = True

NY = 30
N = 10

# plot

if curv_dose:
    # should use same N to emphasise symmetry
    doses = np.linspace(0, 1, N+1)[1:]
    curvatures = np.linspace(0, 10, N+1)[1:]
    rf = 10**(-3)

    x, y, z, conf_str = changing_fcide_dose_curve(doses, curvatures, rf, NY)
    labels = dict(x="Dose", y="Curvature")
    fcide_grid(x, y, z, conf_str, labels)


if curv_asymp:
    # can use different N
    asymps = np.linspace(0, 1, N+1)[1:]
    curvatures = np.linspace(0, 9.6, N+1)[1:]
    rf = 10**(-3)

    x, y, z, conf_str = changing_fcide_curve_asymp(curvatures, asymps, rf, NY)
    labels = dict(x="Curvature", y="Asymptote")
    fcide_grid(x, y, z, conf_str, labels)