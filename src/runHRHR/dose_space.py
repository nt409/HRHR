from utils.functions import RunModel
from utils.plotting import dose_space_contour
from .config import ConfigGridRun

dose_grid_plot = True

ConfigGridRun.config_string = ConfigGridRun.config_string.replace("grid", "dose_space")
output = RunModel().master_loop_dose_space_coordinate(ConfigGridRun)

if dose_grid_plot:
    to_plot = 'LTY'
    fig = dose_space_contour(output, to_plot)
    fig.show()