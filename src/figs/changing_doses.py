import numpy as np

from utils.functions import get_SR_by_doses
from utils.plotting.figures import SR_by_dose_plot


SR_by_dose = True

# setup

doses = np.linspace(0, 1, 10)
freqs = [10**(-5), 0.01, 0.02, 0.05, 0.1]


# plot

if SR_by_dose:
    outputs, conf_str = get_SR_by_doses(doses, freqs)
    SR_by_dose_plot(outputs, conf_str)

