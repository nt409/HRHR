from utils.functions import RunModel
from utils.plotting import yield_by_year, res_freqs_single_t_plot, \
    single_year_plot, yield_res_freqs_plot, plot_frequencies, plot_frequencies_over_time
from .config import ConfigSingleRun

yield_single = False
res_freqs_single = False
yield_res_freqs = False
single_year = False
freq_plot = False
freq_time_plot = True

output = RunModel().master_loop_one_tactic(ConfigSingleRun)
conf_str = ConfigSingleRun.config_string
conf_str = conf_str.replace("saved_runs", "figures")
conf_str = conf_str.replace("pickle", "png")

if yield_single:
    fig = yield_by_year(output)
    fig.show()
    filename = conf_str.replace("/single/", "/single/yield_by_year/")
    fig.write_image(filename)

if yield_res_freqs:
    fig = yield_res_freqs_plot(output)
    fig.show()
    filename = conf_str.replace("/single/", "/single/yield_rf/")
    fig.write_image(filename)

if res_freqs_single:
    fig = res_freqs_single_t_plot(output)
    fig.show()
    filename = conf_str.replace("/single/", "/single/res_freqs/")
    fig.write_image(filename)

if single_year:
    indices_to_plot = list(range(16))
    fig = single_year_plot(output, indices_to_plot)
    fig.show()
    filename = conf_str.replace("/single/", "/single/within_season/")
    fig.write_image(filename)

if freq_plot:
    fig = plot_frequencies(output)
    fig.show()

if freq_time_plot:
    fig = plot_frequencies_over_time(output)
    fig.show()