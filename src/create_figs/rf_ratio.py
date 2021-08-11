from utils.functions import RunRfRatio, RunGrid, compare_me_vs_hobb_by_ratio
from utils.config import ConfigGridRun

from plotting.figures import dose_grid_heatmap_with_contours, outcomes_by_ratio

# which plot

rf_ratio = False
compare_contours = True


# plot

if rf_ratio:
    grid = RunGrid().grid_of_tactics(ConfigGridRun)

    contours = RunRfRatio(grid).get_contours()

    conf_str = ConfigGridRun.config_string_img
    dose_grid_heatmap_with_contours(grid, ConfigGridRun, contours, conf_str)


if compare_contours:
    ratio_list = [10**(i) for i in [-4, -3, -2, -1, 0, 1, 2, 3, 4]]

    output = compare_me_vs_hobb_by_ratio(ConfigGridRun, 10**(-5), ratio_list)
    print(output)
    
    conf_str = ConfigGridRun.config_string_img
    outcomes_by_ratio(output, conf_str)

    # dose_grid_heatmap_with_contours(grid, ConfigGridRun, contours, conf_str)


