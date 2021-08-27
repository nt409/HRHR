from model.simulator import RunSingleTactic, RunGrid

from model.config_classes import SingleConfig, GridConfig

from plotting.figures import DiseaseProgressCurvesAll, DoseSpaceScenariosPlot, \
    DosesScatterPlot, YieldAndRfPlot, ParamScanPlotMeVsHobb, \
    ParamScanPlotHighLowDose, CombinedModelPlot

from param_scan.fns.config import config_rand
from param_scan.fns.post_process import PostProcess



# which plots

model_output_overview = False
rf_yield = False
model_output_combined = False
dose_space = True
doses_scatter = False
param_scan_hobb_vs_me = False
param_scan_high_low_dose = False



def get_param_data(par_str):
    PP = PostProcess(par_str)
    PP.get_maximum_along_contour_df()
    data = PP.max_along_contour_df
    return data


# plot

if model_output_overview:
    config_sing = SingleConfig(1, 2*10**(-1), 5*10**(-2), 1, 1, 0.5, 0.5)
    config_sing.load_saved = False
    output = RunSingleTactic().run(config_sing)
    DiseaseProgressCurvesAll(output, config_sing.config_string_img)


if rf_yield:
    config_sing = SingleConfig(10, 10**(-3), 10**(-5), 1, 1, 0.5, 0.5)
    # config_sing.load_saved = False
    run_s = RunSingleTactic()
    run_s.yield_stopper = 0
    output = run_s.run(config_sing)
    YieldAndRfPlot(output, config_sing.config_string_img)


if model_output_combined:
    config_sing = SingleConfig(10, 10**(-3), 10**(-6), 1, 1, 0.5, 0.5)
    config_sing.load_saved = False
    run_s = RunSingleTactic()
    run_s.yield_stopper = 0
    output = run_s.run(config_sing)
    CombinedModelPlot(output, config_sing.config_string_img)


if dose_space:
    # conf_grid = GridConfig(30, 10**(-7), 10**(-3), 51)

    conf_grid = GridConfig(30, 10**(-7), 10**(-3), 6)
    conf_grid.load_saved = False
    output = RunGrid().run(conf_grid)
    
    DoseSpaceScenariosPlot(output, conf_grid.config_string_img)


if doses_scatter:
    conf_grid = GridConfig(30, 10**(-7), 10**(-3), 51)    
    output = RunGrid().run(conf_grid)
    DosesScatterPlot(output, conf_grid.config_string_img)


if param_scan_hobb_vs_me:
    par_str = config_rand['par_str']
    data = get_param_data(par_str)
    ParamScanPlotMeVsHobb(data, f"{par_str}.png")


if param_scan_high_low_dose:
    par_str = config_rand['par_str']
    data = get_param_data(par_str)
    ParamScanPlotHighLowDose(data, f"{par_str}.png")
