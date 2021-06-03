from utils.functions import RunSingleTactic, RunGrid
from utils.plotting import DiseaseProgressCurvesAll, DoseSpaceScenariosPlot, \
    DosesScatterPlot

from .config_classes import SingleConfig, GridConfig

# which plots

model_output_overview = False
dose_space = False
doses_scatter = True





# plot

if model_output_overview:
    ConfigSingleRun = SingleConfig(1, 2*10**(-1), 5*10**(-2), 1, 1, 0.5, 0.5)

    output = RunSingleTactic().run_single_tactic(ConfigSingleRun)

    DiseaseProgressCurvesAll(output, ConfigSingleRun.config_string_img)

if dose_space:
    ConfigGridRun = GridConfig(30, 10**(-7), 10**(-3), 51)
    output = RunGrid().grid_of_tactics(ConfigGridRun)
    DoseSpaceScenariosPlot(output, ConfigGridRun.config_string_img)

if doses_scatter:
    ConfigGridRun = GridConfig(30, 10**(-7), 10**(-3), 51)
    output = RunGrid().grid_of_tactics(ConfigGridRun)
    DosesScatterPlot(output, ConfigGridRun.config_string_img)

