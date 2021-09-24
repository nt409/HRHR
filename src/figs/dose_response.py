import numpy as np
from math import exp


from model.params import PARAMS
from plotting.paper_figs import DoseResponse
from model.utils import Fungicide


fcd = Fungicide(1, 9.6, PARAMS.delta_1)
fc = Fungicide(0.8, 7, PARAMS.delta_1)

x = np.linspace(0, 1, 100)
xs = [x, x]

ys = [
    [fc.effect(xx) for xx in x],
    [fcd.effect(xx) for xx in x],
    # [1 - fcd.effect(xx) for xx in x],
    # [1 - fc.effect(xx) for xx in x]
    ]


names = ["Alternative weaker fungicide", "Pyraclostrobin"]
cols = ["turquoise", "black"]
dashes = ["solid", "dash"]

dose_resp_data = dict(x=xs, y=ys, name=names, cols=cols, dashes=dashes)


time = np.linspace(1212,2066,1+2066-1212)
concs1 = [exp( - PARAMS.delta_1 * (xx - 1456)) if xx>1456 else 0 for xx in time]
concs2 = [exp( - PARAMS.delta_1 * (xx - 1700)) if xx>1700 else 0 for xx in time]
concs = np.array(concs1) + np.array(concs2)

concs_trcs = [concs, concs]

t1_ind = np.where(time==1456)[0][0]
t2_ind = np.where(time==1700)[0][0]
# print(tt)
# exit()
time = list(time)
time[t1_ind] = None
time[t2_ind] = None
xs = [time, time]

effects = [
    [fc.effect(xx) for xx in concs],
    [fcd.effect(xx) for xx in concs],
    ]

effect_vs_time_data = dict(x=xs, y=effects, concs=concs_trcs, name=names, cols=cols, dashes=dashes)

filename = "../outputs/figures/paper_figs/dose_response.png"

DoseResponse(dose_resp_data, effect_vs_time_data, filename)