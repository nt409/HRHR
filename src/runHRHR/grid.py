from utils.functions import RunModel
from utils.plotting import dose_grid_heatmap, dose_grid_heatmap_with_log_ratio
from .config import ConfigGridRun

# which plots

dose_grid_plot = True
dose_grid_plot_with_log_ratio = True


# run

output = RunModel().master_loop_grid_of_tactics(ConfigGridRun)
conf_str = ConfigGridRun.config_string_img


# plot

if dose_grid_plot:
    to_plot = 'FY'
    dose_grid_heatmap(output, ConfigGridRun, to_plot, conf_str)


if dose_grid_plot_with_log_ratio:
    to_plot = 'LTY'
    dose_grid_heatmap_with_log_ratio(output,
                                ConfigGridRun,
                                to_plot,
                                conf_str)
