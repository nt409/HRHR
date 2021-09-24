import numpy as np

from plotting.paper_figs import DeltaCurve

x = np.linspace(0,1,100)
y = [xx*(1-xx) for xx in x]
y2 = [2*xx*(1-xx) for xx in x]
y3 = [(1-xx**2) for xx in x]

dashes = ["solid", "dash", "dash"]

names = ["Single", "Double (sexual)", "Double (asexual)"]

data = dict(x=[x, x, x], y=[y, y2, y3], dash=dashes, name=names)
filename = "../outputs/figures/paper_figs/delta_curve.png"

DeltaCurve(data, filename)