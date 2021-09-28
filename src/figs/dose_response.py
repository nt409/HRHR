import numpy as np
from math import exp


from model.params import PARAMS
from plotting.paper_figs import DoseResponse
from model.utils import Fungicide


fc = Fungicide(0.8, 5, PARAMS.delta_1)
fcd = Fungicide(1, 9.6, PARAMS.delta_1)

x = np.linspace(0, 1, 100)
xs = [x, x]

ys = [
    [fc.effect(xx) for xx in x],
    [fcd.effect(xx) for xx in x],
    ]


names = ["Alternative weaker fungicide", "Pyraclostrobin"]
cols = ["turquoise", "black"]
dashes = ["solid", "dash"]

dose_resp_data = dict(x=xs, y=ys, name=names, cols=cols, dashes=dashes)

T_end = 2900
time = np.linspace(1212,T_end,1+T_end-1212)
concs1D = [exp( - PARAMS.delta_1 * (xx - 1456)) if xx>1456 else 0 for xx in time]
concs2D = [exp( - PARAMS.delta_1 * (xx - 1700)) if xx>1700 else 0 for xx in time]
concsD = np.array(concs1D) + np.array(concs2D)

concs1A = [exp( - 8e-3 * (xx - 1456)) if xx>1456 else 0 for xx in time]
concs2A = [exp( - 8e-3 * (xx - 1700)) if xx>1700 else 0 for xx in time]
concsA = np.array(concs1A) + np.array(concs2A)

concs_trcs = [concsA, concsD]

t1_ind = np.where(time==1456)[0][0]
t2_ind = np.where(time==1700)[0][0]

time = list(time)
time[t1_ind] = None
time[t2_ind] = None
xs = [time, time]

effects = [
    [fc.effect(xx) for xx in concsA],
    [fcd.effect(xx) for xx in concsD],
    ]

effect_vs_time_data = dict(x=xs, y=effects, concs=concs_trcs, name=names, cols=cols, dashes=dashes)

filename = "../outputs/figures/paper_figs/dose_response.png"

DoseResponse(dose_resp_data, effect_vs_time_data, filename)