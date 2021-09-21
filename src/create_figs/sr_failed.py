"""
Plot failed SR run
"""

from plotting.paper_figs import DoseSpaceScenarioDouble, DoseSpaceScenarioSingle
from model.simulator import RunGrid
from param_scan.fns.config import config_rand
from param_scan.fns.post_process import PostProcess


run_old = False
run_pert = True

PP = PostProcess(config_rand['par_str'])

if run_old:
    PP.re_run_grid(NDoses=201, run_indices=[117])



if run_pert:
    df_test = PP.get_params_for_specific_runs([117])

    pars = df_test.iloc[0,:]

    this_run_ind = int(df_test.iloc[0,:].run)

    rp = PP._get_RPs(pars, 201)
    
    grid_default = RunGrid(rp.fung_parms).run(rp.grid_conf)

    print(f"\nRe-running run: {this_run_ind} \n")

    pars['SR'] = pars['SR']*1.1
    pars['RS'] = pars['RS']*0.9
    pars['RR'] = pars['RR']*1.1

    rp = PP._get_RPs(pars, 201)
    
    grid_pert = RunGrid(rp.fung_parms).run(rp.grid_conf)

    conf_str = rp.grid_conf.config_string_img

    # plot output
    conf_str = conf_str.replace("param_scan/", f"paper_figs/run={this_run_ind}_pert_")
    DoseSpaceScenarioDouble(grid_default, grid_pert, conf_str)


