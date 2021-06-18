from utils.functions import RunDoseSpace, RunGrid, RunRadial

from utils.plotting import dose_space_contour, eq_RFB_contours, radial, \
    dose_sum_LR, dose_sum_hobb_vs_me, first_year_yield, MS_RFB_scatter_plot

from .config import ConfigGridRun

import copy

# which plot

dose_grid_plot = False
dose_sum_hobb_vs_me_plot = False
dose_sum_plot = False
radial_plot = False
first_year_yield_plot = False
RFB_contours = False
MS_RFB_scatter = True


# plot

if dose_grid_plot:
    ConfDGP = copy.copy(ConfigGridRun)

    output = RunDoseSpace().run_loop(ConfDGP)

    to_plot = 'LTY'
    conf_str = ConfDGP.config_string_img
    dose_space_contour(output, to_plot, conf_str)



if dose_sum_hobb_vs_me_plot:
    ConfDS_HM = copy.copy(ConfigGridRun)

    to_plot = 'FY'
    
    ConfDS_HM.n_doses = 20
    ConfDS_HM.res_props = dict(
            f1 = 10**(-3),
            f2 = 10**(-7)
            )

    ConfDS_HM.add_string()
    conf_str = ConfDS_HM.config_string_img

    output = RunGrid().grid_of_tactics(ConfDS_HM)


    dose_sum_hobb_vs_me(output,
                ConfDS_HM,
                to_plot,
                conf_str)


if dose_sum_plot:
    ConfDSP = copy.copy(ConfigGridRun)

    to_plot = 'FY'
    
    ConfDSP.n_doses = 9
    ConfDSP.add_string()
    conf_str = ConfDSP.config_string_img

    output = RunGrid().grid_of_tactics(ConfDSP)

    dose_sum_LR(output,
                ConfDSP,
                to_plot,
                conf_str)


if radial_plot:
    ConfRad = copy.copy(ConfigGridRun)

    ConfRad.n_angles = 5
    ConfRad.n_radii = 10

    ext_string = f"_Nr={ConfRad.n_radii},Na={ConfRad.n_angles}"
    ConfRad.add_string(extra_detail=ext_string)

    radial_output = RunRadial().master_loop_radial(ConfRad)
    
    ConfRad.n_doses = 10

    grid_output = RunGrid().grid_of_tactics(ConfRad)

    radial(radial_output, grid_output, ConfRad)



if first_year_yield_plot:
    ConfFY = copy.copy(ConfigGridRun)

    ConfFY.res_props = dict(
        f1 = 10**(-5),
        f2 = 10**(-3),
        )
    ConfFY.n_doses = 25
    ConfFY.add_string()

    output = RunGrid().grid_of_tactics(ConfFY)

    first_year_yield(output, ConfFY)



if RFB_contours:
    ConfRFB = copy.copy(ConfigGridRun)
    
    ConfRFB.res_props = dict(
        f1 = 10**(-5),
        f2 = 10**(-3),
        )
    ConfRFB.n_doses = 25
    ConfRFB.add_string()

    output = RunGrid().grid_of_tactics(ConfRFB)

    eq_RFB_contours(output, ConfRFB)




if MS_RFB_scatter:
    ConfMS_RFB = copy.copy(ConfigGridRun)

    ConfMS_RFB.res_props = dict(
        f1 = 10**(-5),
        f2 = 10**(-5),
        )
    
    ConfMS_RFB.n_doses = 15
    ConfMS_RFB.add_string()

    output = RunGrid().grid_of_tactics(ConfMS_RFB)

    MS_RFB_scatter_plot(output,
                ConfMS_RFB)
