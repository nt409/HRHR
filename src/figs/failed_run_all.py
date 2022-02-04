"""
Plot failed SR run
"""

import pickle

from plotting.paper_figs import DoseSpace6
from model.simulator import RunGrid
from param_scan.fns.config import config_rand
from param_scan.fns.post_process import PostProcess


def run_all_failed(runs):
    PP = PostProcess(config_rand['par_str'])

    grid_outputs = []

    for run_index in runs:

        filename = f"./param_scan/outputs/failed/grid_output_run={run_index}.pickle"
        with open(filename, 'rb') as f:
            data = pickle.load(f)

        grid_outputs.append(data)

    df_test = PP.get_params_for_specific_runs([run_index])
    pars = df_test.iloc[0, :]
    this_run_ind = int(df_test.iloc[0, :].run)
    rp = PP.get_RPs(pars, 201)
    conf_str = rp.grid_conf.config_string_img

    # plot output
    conf_str = conf_str.replace(
        "param_scan/", f"paper_figs/failed_run_all={this_run_ind}_")
    filename = conf_str[:100] + ".png"
    DoseSpace6(grid_outputs, filename)


if __name__ == "__main__":

    runs = [
        91,
        183,
        227,
        241,
        295,
        478,
    ]

    run_all_failed(runs)
