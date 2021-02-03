from utils.functions import RunModel
from utils.plotting import yield_by_year, res_freqs_single_t_plot, \
    single_year_plot
from .config import ConfigSingleRun

yield_single = True
res_freqs_single = True
single_year = True

output = RunModel().master_loop_one_tactic(ConfigSingleRun)

if yield_single:
    fig = yield_by_year(output)
    fig.show()

if res_freqs_single:
    fig = res_freqs_single_t_plot(output)
    fig.show()

if single_year:
    indices_to_plot = list(range(16))
    fig = single_year_plot(output, indices_to_plot)
    fig.show()