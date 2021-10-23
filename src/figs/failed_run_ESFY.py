"""
Plot failed SR run
"""

import pickle

from plotting.paper_figs import DoseSpaceScenarioSingle
from model.simulator import RunGrid
from param_scan.fns.config import config_rand
from param_scan.fns.post_process import PostProcess



def run_single(run_index):
    df_test = PP.get_params_for_specific_runs([run_index])
    pars = df_test.iloc[0,:]
    this_run_ind = int(df_test.iloc[0,:].run)
    rp = PP._get_RPs(pars, 201)
    conf_str = rp.grid_conf.config_string_img
    
    
    grid_default = RunGrid(rp.fung_parms).run(rp.grid_conf)
    # filename = f"./param_scan/outputs/failed/grid_output_run={this_run_ind}.pickle"
    # with open(filename, 'rb') as f:
        # grid_default = pickle.load(f)


    # plot output
    conf_str = conf_str.replace("param_scan/", f"paper_figs/failed_run_ESFY={this_run_ind}_")
    filename = conf_str[:100] + ".png"
    DoseSpaceScenarioSingle(grid_default, filename)


if __name__=="__main__":

    PP = PostProcess(config_rand['par_str'])

    run_single(265)