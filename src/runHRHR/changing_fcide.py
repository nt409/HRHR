import numpy as np

from utils.functions import changing_fcide_dose_curve
from utils.plotting import fcide_grid

curv_dose = True

# setup

N = 2
doses = np.linspace(0, 1, N+1)[1:]
curvatures = np.linspace(0, 10, N+1)[1:]
rf = 10**(-3)




# plot

if curv_dose:
    x, y, z, conf_str = changing_fcide_dose_curve(doses, curvatures, rf)
    fcide_grid(x, y, z, conf_str)
    
