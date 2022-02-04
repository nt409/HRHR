import numpy as np
from math import exp


from model.params import PARAMS
from plotting.paper_figs import DoseResponse
from model.utils import Fungicide


fc1 = Fungicide(0.8, 5, PARAMS.delta_1)
fc2 = Fungicide(0.8, 10, PARAMS.delta_1)

x = np.linspace(0, 1, 100)
xs = [x, x]

ys = [
    [fc1.effect(xx) for xx in x],
    [fc2.effect(xx) for xx in x],
]


names = ["Low efficacy (curvature 5)", "High efficacy (curvature 10)"]
cols = ["turquoise", "black"]
dashes = ["solid", "dash"]

dose_resp_data = dict(x=xs, y=ys, name=names, cols=cols, dashes=dashes)

T_end = 2900
time = np.linspace(1212, T_end, 1+T_end-1212)

concs1f = [exp(- PARAMS.delta_1 * (xx - 1456))
           if xx > 1456 else 0 for xx in time]
concs1s = [exp(- PARAMS.delta_1 * (xx - 1700))
           if xx > 1700 else 0 for xx in time]
concs1 = np.array(concs1f) + np.array(concs1s)

concs2f = [0.5*exp(- PARAMS.delta_1 * (xx - 1456))
           if xx > 1456 else 0 for xx in time]
concs2s = [0.5*exp(- PARAMS.delta_1 * (xx - 1700))
           if xx > 1700 else 0 for xx in time]
concs2 = np.array(concs2f) + np.array(concs2s)

concs_trcs = [concs1, concs2]

t1_ind = np.where(time == 1456)[0][0]
t2_ind = np.where(time == 1700)[0][0]

time = list(time)
time[t1_ind] = None
time[t2_ind] = None
xs = [time, time]

effects = [
    [fc1.effect(xx) for xx in concs1],
    [fc2.effect(xx) for xx in concs2],
]

effect_vs_time_data = dict(
    x=xs, y=effects, concs=concs_trcs, name=names, cols=cols, dashes=dashes)

filename = "../outputs/figures/paper_figs/dose_response_app.png"

DoseResponse(dose_resp_data, effect_vs_time_data, filename, True)
