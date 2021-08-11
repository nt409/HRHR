from utils.functions import RunGrid
from utils.config import ConfigGridRun

from plotting.figures import dose_grid_heatmap, dose_grid_RA_heatmap

# which plots

dose_grid_plot = False
res_array_plot = True



# plot

if dose_grid_plot:
    to_plot = 'econ'
    output = RunGrid().grid_of_tactics(ConfigGridRun)

    conf_str = ConfigGridRun.config_string_img
    dose_grid_heatmap(output, ConfigGridRun, to_plot, conf_str)


if res_array_plot:
    output = RunGrid().grid_of_tactics(ConfigGridRun)

    conf_str = ConfigGridRun.config_string_img
    for i in range(15):
        dose_grid_RA_heatmap(output, ConfigGridRun, conf_str, i)

