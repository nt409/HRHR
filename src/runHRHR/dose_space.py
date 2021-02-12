from utils.functions import RunModel
from utils.plotting import dose_space_contour, radial
from .config import ConfigGridRun

# which plot

dose_grid_plot = False
radial_plot = True

# plot

if dose_grid_plot:
    output = RunModel().master_loop_dose_space_coordinate(ConfigGridRun)

    to_plot = 'LTY'
    conf_str = ConfigGridRun.config_string_img
    dose_space_contour(output, to_plot, conf_str)
    

if radial_plot:
    ConfigGridRun.n_angles = 5
    ConfigGridRun.n_radii = 7

    radial_output = RunModel().master_loop_radial(ConfigGridRun)
    
    ConfigGridRun.n_doses = 10
    ConfigGridRun.add_string()
    grid_output = RunModel().master_loop_grid_of_tactics(ConfigGridRun)

    radial(radial_output, grid_output, ConfigGridRun)