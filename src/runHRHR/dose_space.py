from utils.functions import RunModel
from utils.plotting import dose_space_contour, radial
from .config import ConfigGridRun

dose_grid_plot = True
radial_plot = True

conf_str = ConfigGridRun.config_string
conf_str = conf_str.replace("saved_runs", "figures")
conf_str = conf_str.replace("pickle", "png")

if dose_grid_plot:
    output = RunModel().master_loop_dose_space_coordinate(ConfigGridRun)

    to_plot = 'LTY'
    fig = dose_space_contour(output, to_plot)
    fig.show()
    filename = conf_str.replace("/grid/", "/dose_space/yield_by_year/")
    fig.write_image(filename)

if radial_plot:
    radial_output = RunModel().master_loop_radial(ConfigGridRun)
    
    ConfigGridRun.n_doses = 10
    ConfigGridRun.add_string()
    grid_output = RunModel().master_loop_grid_of_tactics(ConfigGridRun)

    fig = radial(radial_output, grid_output, ConfigGridRun)
    fig.show()
    filename = conf_str.replace("/grid/", "/dose_space/radial/")
    fig.write_image(filename)