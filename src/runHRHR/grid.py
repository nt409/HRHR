from utils.functions import RunModel
from utils.plotting import dose_grid_heatmap
from .config import ConfigGridRun

dose_grid_plot = True


output = RunModel().master_loop_grid_of_tactics(ConfigGridRun)

if dose_grid_plot:
    to_plot = 'FY'
    fig = dose_grid_heatmap(output, ConfigGridRun, to_plot)
    fig.show()