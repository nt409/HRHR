from utils.functions import RunModel
from utils.plotting import yield_by_year, res_freqs_single_t_plot, \
    single_year_plot, yield_res_freqs_plot, plot_frequencies, plot_frequencies_over_time
from .config import ConfigSingleRun

# which plots

yield_single = False
res_freqs_single = False
yield_res_freqs = True
single_year = False
freq_plot = False
freq_time_plot = False

# run

output = RunModel().master_loop_one_tactic(ConfigSingleRun)
conf_str = ConfigSingleRun.config_string_img


# plots

if yield_single:
    yield_by_year(output, conf_str)


if yield_res_freqs:
    yield_res_freqs_plot(output, conf_str)
    

if res_freqs_single:
    res_freqs_single_t_plot(output, conf_str)


if single_year:
    indices_to_plot = list(range(16))
    single_year_plot(output, indices_to_plot, conf_str)


if freq_plot:
    plot_frequencies(output, conf_str)
    

if freq_time_plot:
    plot_frequencies_over_time(output, conf_str)

