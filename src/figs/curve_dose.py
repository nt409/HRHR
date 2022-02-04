"""NB now using dose response app!"""

import numpy as np
from math import exp


from model.params import PARAMS
from plotting.paper_figs import CurveDose
from model.utils import Fungicide


fcl = Fungicide(0.8, 5, PARAMS.delta_1)
fch = Fungicide(0.8, 10, PARAMS.delta_1)

xA = np.linspace(0, 2, 100)
xB = np.linspace(0, 1, 100)
xs = [xA, xB]

ys = [
    [fcl.effect(xx) for xx in xA],
    [fch.effect(xx) for xx in xB],
]


names = ["Curvature 5", "Curvature 10"]
cols = ["turquoise", "black"]
dashes = ["solid", "dash"]

dose_resp_data = dict(x=xs, y=ys, name=names, cols=cols, dashes=dashes)


time = np.linspace(0, 600, 601)
concsA1 = [exp(- PARAMS.delta_1 * xx) for xx in time]
concsA2 = [exp(- PARAMS.delta_1 * (xx - 244))
           if xx > 244 else 0 for xx in time]
concsA = np.array(concsA1) + np.array(concsA2)

concsB1 = [0.5 * exp(- PARAMS.delta_1 * xx) for xx in time]
concsB2 = [0.5 * exp(- PARAMS.delta_1 * (xx - 244))
           if xx > 244 else 0 for xx in time]
concsB = np.array(concsB1) + np.array(concsB2)

xs = [time, time]

effects = [
    [fcl.effect(xx) for xx in concsA],
    [fch.effect(xx) for xx in concsB],
]

effect_vs_time_data = dict(
    x=xs, y=effects, name=names, cols=cols, dashes=dashes)

filename = "../outputs/figures/paper_figs/curve_dose.png"

CurveDose(dose_resp_data, effect_vs_time_data, filename)
