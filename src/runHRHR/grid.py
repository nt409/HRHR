from utils.functions import RunModel
from utils.plotting import dose_grid_heatmap
from .config import ConfigGridRun

# which plots

dose_grid_plot = False



# plot

if dose_grid_plot:
    to_plot = 'FY'
    output = RunModel().grid_of_tactics(ConfigGridRun)

    conf_str = ConfigGridRun.config_string_img
    dose_grid_heatmap(output, ConfigGridRun, to_plot, conf_str)

