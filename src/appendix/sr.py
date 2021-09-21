from math import exp, log
from scipy.integrate import simps
import numpy as np

from model.params import PARAMS
import plotly.graph_objects as go


x = log(2) / 9.6
print(x)
exit()

def fung(t, dose):
    out = 1 - exp( - PARAMS.theta_1 * dose * exp(-(PARAMS.delta_1)*t))
    # out = 1 - dose
    return out

rho = 1 * PARAMS.beta

def get_ysr(dose, ts):
    return [exp(rho*fung(t, dose)) for t in ts]

def get_yss(dose, ts):
    return [exp(rho*(fung(t, dose)**2)) for t in ts]


ts = list(range(150))
doses = np.linspace(0,1,10)

y = []

for dose in doses:
    ysr = get_ysr(dose, ts)
    yss = get_yss(dose, ts)

    int_sr = simps(ysr, ts)
    int_ss = simps(yss, ts)

    R = int_sr/int_ss
    y.append(R)


trc = go.Scatter(x=doses, y=y)
fig = go.Figure(data=trc)
fig.show()