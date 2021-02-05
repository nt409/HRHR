from utils.functions import RunModel
from utils.plotting import dose_grid_heatmap, dose_grid_heatmap_with_log_ratio
from .config import ConfigGridRun

dose_grid_plot = True
dose_grid_plot_with_log_ratio = True

output = RunModel().master_loop_grid_of_tactics(ConfigGridRun)
conf_str = ConfigGridRun.config_string
conf_str = conf_str.replace("saved_runs", "figures")
conf_str = conf_str.replace("pickle", "png")

if dose_grid_plot:
    to_plot = 'FY'
    fig = dose_grid_heatmap(output, ConfigGridRun, to_plot)
    fig.show()
    filename = conf_str.replace("/grid/", "/grid/dose_grid/")
    fig.write_image(filename)

if dose_grid_plot_with_log_ratio:
    to_plot = 'LTY'
    fig = dose_grid_heatmap_with_log_ratio(output, ConfigGridRun, to_plot)
    fig.show()
    filename = conf_str.replace("/grid/", "/grid/dose_grid_LR/")
    fig.write_image(filename)