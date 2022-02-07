import scipy
import sys

from param_scan.run.r4_re_run_failures import check_outcome_RFB

run_attrs = [
    dict(run=10,  DS_lim=[0.5, 1.1]),  # Y
    dict(run=116, DS_lim=[0.2, 1.0]),  # Y
    dict(run=295, DS_lim=[0.67, 0.7]),  # Y DS_lim=[0.6,0.75]

    dict(run=91,  DS_lim=[0.4, 0.5]),  # N
    dict(run=183, DS_lim=[0.3, 0.45]),  # N
    dict(run=227, DS_lim=[0.6, 0.7]),  # M?
    dict(run=241, DS_lim=[0.3, 0.6]),  # N

    dict(run=478, DS_lim=[0.2, 0.6]),  # N

    # ESFY outperformed ERFB
    # dict(run=265, DS_lim=[0.44, 0.46]),  # N
]


if __name__ == "__main__":

    if len(sys.argv) != 2:
        print(len(sys.argv))
        raise Exception("Supply one argument: the run index")

    index = int(sys.argv[1])
    run_attrs_use = [run_attrs[index]]

    N_grid_doses = 101
    N_cont_doses = 500

    run_attrs_use[0]["NDoses"] = N_grid_doses
    run_attrs_use[0]["N_cont_doses"] = N_cont_doses

    check_outcome_RFB(run_attrs_use)
    # check_outcome_SFY(run_attrs_use)
