from utils.functions import RunGrid
from utils.plotting import dose_grid_heatmap
from .config import ConfigGridRun

# which plots

dose_grid_plot = True



# plot

if dose_grid_plot:
    to_plot = 'econ'
    output = RunGrid().grid_of_tactics(ConfigGridRun)

    conf_str = ConfigGridRun.config_string_img
    dose_grid_heatmap(output, ConfigGridRun, to_plot, conf_str)

