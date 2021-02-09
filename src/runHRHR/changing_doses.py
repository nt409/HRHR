import itertools
import numpy as np

from utils.functions import RunModel
from utils.plotting import SR_by_dose_plot
from .config_classes import SingleConfig

SR_by_dose = True

doses = np.linspace(0, 1, 10)
freqs = [10**(-5), 0.01, 0.02, 0.05, 0.1, 0.2]




# run

outputs = {}
for dose, rf in itertools.product(doses, freqs):
    ConfigSingleRun = SingleConfig(1, rf, rf, dose, dose, dose, dose)
    output = RunModel().master_loop_one_tactic(ConfigSingleRun)
    outputs[f"dose={dose},rf={rf}"] = output


conf_str = ConfigSingleRun.config_string
str_freqs = [str(round(f,2)) for f in freqs]
str_doses = [str(round(d,2)) for d in doses]

middle_string = ("=" + ",_".join(str_freqs) +
             "_doses=" + ",_".join(str_doses))
middle_string = middle_string.replace(".", ",")

conf_str = ("=".join(conf_str.split("=")[0:2]) + 
        middle_string + conf_str.split("=")[-1])

conf_str = conf_str.replace("saved_runs", "figures")
conf_str = conf_str.replace("pickle", "png")


# plot

if SR_by_dose:
    fig = SR_by_dose_plot(outputs)
    fig.show()
    filename = conf_str.replace("/single/", "/changing_dose/equal_ratio/")
    fig.write_image(filename)
