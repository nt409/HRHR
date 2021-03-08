from utils.functions import RunSingleTactic
from utils.plotting import yield_by_year, res_freqs_single_t_plot, \
    single_year_plot, yield_res_freqs_plot, plot_frequencies, plot_frequencies_over_time
from .config import ConfigSingleRun
from .config_classes import SingleConfig


# which plots

yield_single = False
res_freqs_single = False
yield_res_freqs = False
single_year = True
freq_bar_plot = False
freq_time_plot = False

# run
bools = [res_freqs_single, freq_bar_plot]

if any(bools):
    output = RunSingleTactic().run_single_tactic(ConfigSingleRun)


# plots
conf_str = ConfigSingleRun.config_string_img

if yield_single:
    rf1 = 10**(-4)
    rf2 = 10**(-2)
    ConfigYRF = SingleConfig(20, rf1, rf2, 1, 1, 1, 1)

    output = RunSingleTactic().run_single_tactic(ConfigYRF)

    yield_by_year(output, conf_str)


if yield_res_freqs:
    rf1 = 10**(-4)
    rf2 = 10**(-2)
    ConfigYRF = SingleConfig(20, rf1, rf2, 1, 1, 1, 1)

    output = RunSingleTactic().run_single_tactic(ConfigYRF)

    yield_res_freqs_plot(output, ConfigYRF.config_string_img)
    

if res_freqs_single:
    res_freqs_single_t_plot(output, conf_str)


if single_year:
    rf1 = 10**(-1)
    rf2 = 2*10**(-1)
    ConfigSY = SingleConfig(20, rf1, rf2, 1, 1, 0.5, 0.5)
    output = RunSingleTactic().run_single_tactic(ConfigSY)

    indices_use = dict(
        # all = list(range(16)),
        sus_rec = [0, 9],
        infected = list(range(5,9)),
        fungs = [14,15],
        )
    
    for key in indices_use:
        indices_to_plot = indices_use[key]
        single_year_plot(output, indices_to_plot, conf_str)


if freq_bar_plot:
    plot_frequencies(output, conf_str)
    

if freq_time_plot:
    rf1 = 10**(-4)
    rf2 = 10**(-2)
    ConfigFTP = SingleConfig(10, rf1, rf2, 1, 1, 1, 1)
    ConfigFTP.sex_prop = 0
    ConfigFTP.add_string()

    output = RunSingleTactic().run_single_tactic(ConfigFTP)

    plot_frequencies_over_time(output, ConfigFTP.config_string_img)

