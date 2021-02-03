import itertools
import numpy as np

from utils.functions import RunModel
from utils.plotting import SR_by_dose_plot
from .config_classes import SingleConfig

SR_by_dose = True

doses = np.linspace(0, 1, 10)
freqs = [10**(-5), 0.01, 0.02, 0.05, 0.1]

# run

outputs = {}
for dose, rf in itertools.product(doses, freqs):
    ConfigSingleRun = SingleConfig(1, rf, rf, dose, dose, dose, dose)
    output = RunModel().master_loop_one_tactic(ConfigSingleRun)
    outputs[f"dose={dose},rf={rf}"] = output


# plot

if SR_by_dose:
    fig = SR_by_dose_plot(outputs)
    fig.show()
