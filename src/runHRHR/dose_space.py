from utils.functions import RunModel
from utils.plotting import dose_space_contour, radial, \
    dose_sum_LR
from .config import ConfigGridRun

# which plot

dose_grid_plot = False
radial_plot = False
dose_sum_plot = True


# plot

if dose_grid_plot:
    output = RunModel().master_loop_dose_space_coordinate(ConfigGridRun)

    to_plot = 'LTY'
    conf_str = ConfigGridRun.config_string_img
    dose_space_contour(output, to_plot, conf_str)



if dose_sum_plot:
    to_plot = 'FY'
    
    ConfigGridRun.n_doses = 13
    ConfigGridRun.add_string()
    conf_str = ConfigGridRun.config_string_img

    output = RunModel().grid_of_tactics(ConfigGridRun)


    dose_sum_LR(output,
                ConfigGridRun,
                to_plot,
                conf_str)


if radial_plot:
    ConfigGridRun.n_angles = 5
    ConfigGridRun.n_radii = 10

    ext_string = f"_Nr={ConfigGridRun.n_radii},Na={ConfigGridRun.n_angles}"
    ConfigGridRun.add_string(extra_detail=ext_string)

    radial_output = RunModel().master_loop_radial(ConfigGridRun)
    
    ConfigGridRun.n_doses = 10

    grid_output = RunModel().grid_of_tactics(ConfigGridRun)

    radial(radial_output, grid_output, ConfigGridRun)