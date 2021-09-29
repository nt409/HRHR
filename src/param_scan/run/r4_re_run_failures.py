"""
Re-run failed cases - try to find out why
"""

import sys

from param_scan.fns.config import config_rand
from param_scan.fns.post_process import PostProcess


def main(run_attrs):
        
    PP = PostProcess(config_rand['par_str'])


    df = PP.processed_df
    print("FAILURES:")
    print(df.loc[df["c_R_maxCont%"]<100, ["max_grid_EL", "c_R_maxContEL", "I_R_best_value"]])

    print("ESFY 'outpeforms' ERFB:")
    print(df.loc[df["c_R_maxContEL"]<df["c_E_maxContEL"], ["max_grid_EL", "c_E_maxContEL", "c_R_maxContEL", "I_R_best_value"]])


    for data in run_attrs:
        NDoses = data["NDoses"]
        N_cont_doses = data["N_cont_doses"]
        DS_lim = data["DS_lim"]
        run = data["run"]        

        # saves result to outpus/re_run/...
        PP.re_run_cont(NDoses=NDoses, N_cont_doses=N_cont_doses, DS_lim=DS_lim, run_indices=[run])
        # PP.re_run_grid(NDoses=NDoses, run_indices=[run])




if __name__=="__main__":
    
    N_dose_def = 51
    N_cont_ds = 100
    
    # N_dose_def = 11
    # N_cont_ds = 5

    run_attrs = [
        dict(run = 10,  DS_lim=[0.5,1.1], NDoses=N_dose_def, N_cont_doses=N_cont_ds),
        dict(run = 91,  DS_lim=[0.3,1.0], NDoses=N_dose_def, N_cont_doses=N_cont_ds),
        dict(run = 116, DS_lim=[0.2,1.0], NDoses=N_dose_def, N_cont_doses=N_cont_ds),
        dict(run = 183, DS_lim=[0.3,0.6], NDoses=N_dose_def, N_cont_doses=N_cont_ds),
        dict(run = 227, DS_lim=[0.5,1.1], NDoses=N_dose_def, N_cont_doses=N_cont_ds),
        dict(run = 241, DS_lim=[0.3,0.6], NDoses=N_dose_def, N_cont_doses=N_cont_ds),
        dict(run = 295, DS_lim=[0.6,1.0], NDoses=N_dose_def, N_cont_doses=N_cont_ds),
        dict(run = 478, DS_lim=[0.2,0.6], NDoses=N_dose_def, N_cont_doses=N_cont_ds),

        # ESFY outperformed ERFB
        dict(run = 265, DS_lim=[0.4,0.8], NDoses=N_dose_def, N_cont_doses=N_cont_ds),
        ]
    
    
    
    if len(sys.argv)!=2:
        raise Exception("Supply one argument: the run index")

    index = int(sys.argv[1])
    
    run_attrs_use = [run_attrs[index]]
    main(run_attrs_use)


