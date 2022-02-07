"""
Re-run failed cases - try to find out why
"""

import sys

# import numpy as np
# from model.strategy_arrays import EqualResFreqBreakdownArray

from model.utils import object_dump

from param_scan.fns.config import config_rand
from param_scan.fns.post_process import PostProcess


def print_it():
    PP = PostProcess(config_rand['par_str'])

    df = PP.processed_df

    print(f"n_runs = {df.shape[0]}")

    print("FAILURES:")
    ERFB_sub_opt = df.loc[df["c_R_maxCont%"] < 100, [
        "max_grid_EL", "c_R_maxContEL", "I_R_best_value", "c_E_maxContEL"]]
    print(ERFB_sub_opt)
    print(f"n failures: {ERFB_sub_opt.shape[0]}")

    print("ESFY 'outperforms' ERFB:")
    ESFY_outperforms = df.loc[df["c_R_maxContEL"] < df["c_E_maxContEL"], [
        "max_grid_EL", "c_E_maxContEL", "c_R_maxContEL", "I_R_best_value"]]
    # print(df.loc[df["c_R_maxContEL"]<df["c_E_maxContEL"]].iloc[0])
    print(ESFY_outperforms)
    print(f"n failures: {ESFY_outperforms.shape[0]}")


def check_outcome_RFB(run_attrs):

    PP = PostProcess(config_rand['par_str'])

    for data in run_attrs:
        N_grid = data["NDoses"]
        run = data["run"]

        # saves result to outputs/re_run/...
        N_cont_doses = data["N_cont_doses"]
        DS_limits = data["DS_lim"]

        PP.re_run_cont_RFB(
            N_grid_doses=N_grid,
            N_cont_doses=N_cont_doses,
            dose_sum_limits=DS_limits,
            run_indices=[run]
        )


def check_outcome_SFY(run_attrs):

    PP = PostProcess(config_rand['par_str'])

    for data in run_attrs:
        N_grid = data["NDoses"]
        run = data["run"]

        # saves result to outputs/re_run/...
        N_cont_doses = data["N_cont_doses"]
        DS_limits = data["DS_lim"]

        PP.re_run_cont_SFY(
            N_grid_doses=N_grid,
            N_cont_doses=N_cont_doses,
            dose_sum_limits=DS_limits,
            run_indices=[run]
        )


def re_run_grid_RFB(run_attrs, plot):

    PP = PostProcess(config_rand['par_str'])

    for data in run_attrs:
        NDoses = data["NDoses"]
        run = data["run"]

        grid_out = PP.re_run_grid(NDoses=NDoses, run_indices=[run], plot=plot)
        filename = f"./param_scan/outputs/failed/grid_output_run={run}.pickle"
        object_dump(filename, grid_out)


if __name__ == "__main__":

    run_attrs = [

        # dict(run=10,  DS_lim=[0.5, 1.1]),  # Y
        # dict(run=116, DS_lim=[0.2, 1.0]),  # Y
        dict(run=295, DS_lim=[0.67, 0.7]),  # Y DS_lim=[0.6,0.75]

        # dict(run=91,  DS_lim=[0.4, 0.5]),  # N
        # dict(run=183, DS_lim=[0.3, 0.45]),  # N
        # dict(run=227, DS_lim=[0.6, 0.7]),  # M?
        # dict(run=241, DS_lim=[0.3, 0.6]),  # N

        # dict(run=478, DS_lim=[0.2, 0.6]),  # N

        # ESFY outperformed ERFB
        # dict(run=265, DS_lim=[0.44, 0.46]),  # N
    ]

    N_grid_doses = 101
    N_cont_doses = 250

    if len(sys.argv) != 2:
        raise Exception("Supply one argument: the run index")

    index = int(sys.argv[1])
    run_attrs_use = [run_attrs[index]]

    run_attrs_use[0]["NDoses"] = N_grid_doses
    run_attrs_use[0]["N_cont_doses"] = N_cont_doses

    # print_it()
    check_outcome_RFB(run_attrs_use)
    # check_outcome_SFY(run_attrs_use)

    # ! to re run failures via cluster
    # re_run_grid_RFB(run_attrs_use, plot=False)
