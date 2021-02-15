from utils.functions import RunModel
from utils.plotting import yield_by_year, res_freqs_single_t_plot, \
    single_year_plot, yield_res_freqs_plot, plot_frequencies, plot_frequencies_over_time
from .config import ConfigSingleRun
from .config_classes import SingleConfig


# which plots

yield_single = True
res_freqs_single = False
yield_res_freqs = False
single_year = False
freq_plot = False
freq_time_plot = False

# run
bools = [yield_single, res_freqs_single, single_year, freq_plot]

if any(bools):
    output = RunModel().run_single_tactic(ConfigSingleRun)


# plots
conf_str = ConfigSingleRun.config_string_img

if yield_single:
    yield_by_year(output, conf_str)


if yield_res_freqs:
    rf1 = 10**(-4)
    rf2 = 10**(-2)
    ConfigYRF = SingleConfig(20, rf1, rf2, 1, 1, 1, 1)

    output = RunModel().run_single_tactic(ConfigYRF)

    yield_res_freqs_plot(output, ConfigYRF.config_string_img)
    

if res_freqs_single:
    res_freqs_single_t_plot(output, conf_str)


if single_year:
    indices_to_plot = list(range(16))
    single_year_plot(output, indices_to_plot, conf_str)


if freq_plot:
    plot_frequencies(output, conf_str)
    

if freq_time_plot:
    rf1 = 10**(-4)
    rf2 = 10**(-2)
    ConfigFTP = SingleConfig(10, rf1, rf2, 1, 1, 1, 1)
    ConfigFTP.sex_prop = 0
    ConfigFTP.add_string()

    output = RunModel().run_single_tactic(ConfigFTP)

    plot_frequencies_over_time(output, ConfigFTP.config_string_img)

