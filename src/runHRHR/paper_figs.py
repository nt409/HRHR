from utils.functions import RunSingleTactic, RunGrid
from utils.plotting import DiseaseProgressCurvesAll, DoseSpaceScenariosPlot, \
    DosesScatterPlot, YieldAndRfPlot, ParamScanPlotMeVsHobb, \
    ParamScanPlotHighLowDose, CombinedModelPlot

from .config_classes import SingleConfig, GridConfig


from param_scan_cluster.config import config_rand
from param_scan_cluster.functions import get_PS_rand_str, PostProcess

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
    ConfigSingleRun = SingleConfig(1, 2*10**(-1), 5*10**(-2), 1, 1, 0.5, 0.5)
    ConfigSingleRun.load_saved = False
    output = RunSingleTactic().run_single_tactic(ConfigSingleRun)
    DiseaseProgressCurvesAll(output, ConfigSingleRun.config_string_img)


if rf_yield:
    ConfigSingleRun = SingleConfig(10, 10**(-3), 10**(-5), 1, 1, 0.5, 0.5)
    # ConfigSingleRun.load_saved = False
    RST = RunSingleTactic()
    RST.yield_stopper = 0
    output = RST.run_single_tactic(ConfigSingleRun)
    YieldAndRfPlot(output, ConfigSingleRun.config_string_img)


if model_output_combined:
    ConfigSingleRun = SingleConfig(10, 10**(-3), 10**(-6), 1, 1, 0.5, 0.5)
    # ConfigSingleRun.load_saved = False
    RST = RunSingleTactic()
    RST.yield_stopper = 0
    output = RST.run_single_tactic(ConfigSingleRun)
    CombinedModelPlot(output, ConfigSingleRun.config_string_img)


if dose_space:
    ConfigGridRun = GridConfig(30, 10**(-7), 10**(-3), 51)
    output = RunGrid().grid_of_tactics(ConfigGridRun)
    DoseSpaceScenariosPlot(output, ConfigGridRun.config_string_img)


if doses_scatter:
    ConfigGridRun = GridConfig(30, 10**(-7), 10**(-3), 51)
    output = RunGrid().grid_of_tactics(ConfigGridRun)
    DosesScatterPlot(output, ConfigGridRun.config_string_img)


if param_scan_hobb_vs_me:
    par_str = get_PS_rand_str(config_rand)
    data = get_param_data(par_str)
    ParamScanPlotMeVsHobb(data, f"{par_str}.png")


if param_scan_high_low_dose:
    par_str = get_PS_rand_str(config_rand)
    data = get_param_data(par_str)
    ParamScanPlotHighLowDose(data, f"{par_str}.png")
