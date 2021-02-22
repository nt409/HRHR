from utils.functions import RunDoseSpace, RunGrid, RunRadial
from utils.plotting import dose_space_contour, radial, \
    dose_sum_LR, dose_sum_hobb_vs_me, first_year_yield
from .config import ConfigGridRun

# which plot

dose_grid_plot = False
dose_sum_hobb_vs_me_plot = False
dose_sum_plot = False
radial_plot = False
first_year_yield_plot = True


# plot

if dose_grid_plot:
    output = RunDoseSpace().run_loop(ConfigGridRun)

    to_plot = 'LTY'
    conf_str = ConfigGridRun.config_string_img
    dose_space_contour(output, to_plot, conf_str)



if dose_sum_hobb_vs_me_plot:
    to_plot = 'FY'
    
    ConfigGridRun.n_doses = 14
    ConfigGridRun.add_string()
    conf_str = ConfigGridRun.config_string_img

    output = RunGrid().grid_of_tactics(ConfigGridRun)


    dose_sum_hobb_vs_me(output,
                ConfigGridRun,
                to_plot,
                conf_str)

if dose_sum_plot:
    to_plot = 'FY'
    
    ConfigGridRun.n_doses = 14
    ConfigGridRun.add_string()
    conf_str = ConfigGridRun.config_string_img

    output = RunGrid().grid_of_tactics(ConfigGridRun)


    dose_sum_LR(output,
                ConfigGridRun,
                to_plot,
                conf_str)


if radial_plot:
    ConfigGridRun.n_angles = 5
    ConfigGridRun.n_radii = 10

    ext_string = f"_Nr={ConfigGridRun.n_radii},Na={ConfigGridRun.n_angles}"
    ConfigGridRun.add_string(extra_detail=ext_string)

    radial_output = RunRadial().master_loop_radial(ConfigGridRun)
    
    ConfigGridRun.n_doses = 10

    grid_output = RunGrid().grid_of_tactics(ConfigGridRun)

    radial(radial_output, grid_output, ConfigGridRun)



if first_year_yield_plot:
    output = RunGrid().grid_of_tactics(ConfigGridRun)

    first_year_yield(output, ConfigGridRun)