"""
Plot failed SR run
"""

import pickle

from plotting.paper_figs import DoseSpaceScenarioDouble
from model.simulator import RunGrid
from param_scan.fns.config import config_rand
from param_scan.fns.post_process import PostProcess


def run_pert(run_index):
    df_test = PP.get_params_for_specific_runs([run_index])

    pars = df_test.iloc[0, :]

    this_run_ind = int(df_test.iloc[0, :].run)

    # rp = PP._get_RPs(pars, NDoses=101)
    # grid_default = RunGrid(rp.fung_parms).run(rp.grid_conf)
    filename = f"./param_scan/outputs/failed/grid_output_run={this_run_ind}.pickle"
    with open(filename, 'rb') as f:
        grid_default = pickle.load(f)

    print(f"\nRe-running run: {this_run_ind} \n")

    pars['SR'] = pars['SR']*1.1
    pars['RS'] = pars['RS']*0.9
    pars['RR'] = pars['RR']*1.1

    rp = PP.get_RPs(pars, NDoses=101)

    grid_pert = RunGrid(rp.fung_parms).run(rp.grid_conf)

    conf_str = rp.grid_conf.config_string_img

    # plot output
    conf_str = conf_str.replace(
        "param_scan/", f"paper_figs/failed_run_pert={this_run_ind}_")
    filename = conf_str[:100] + ".png"
    DoseSpaceScenarioDouble(grid_default, grid_pert, filename)


if __name__ == "__main__":

    run = 91

    PP = PostProcess(config_rand['par_str'])

    run_pert(run)
