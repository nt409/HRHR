"""
Re-run failed cases - try to find out why
"""

import sys
import numpy as np
# from model.strategy_arrays import EqualResFreqBreakdownArray
from model.utils import object_dump

from param_scan.fns.config import config_rand
from param_scan.fns.post_process import PostProcess





def main(run_attrs, plot, check_if_failure, print_it):
        
    PP = PostProcess(config_rand['par_str'])

    if print_it:
        df = PP.processed_df
        print("FAILURES:")
        print(df.loc[df["c_R_maxCont%"]<100, ["max_grid_EL", "c_R_maxContEL", "I_R_best_value"]])

        print("ESFY 'outperforms' ERFB:")
        print(df.loc[df["c_R_maxContEL"]<df["c_E_maxContEL"], ["max_grid_EL", "c_E_maxContEL", "c_R_maxContEL", "I_R_best_value"]])

    for data in run_attrs:
        NDoses = data["NDoses"]
        run = data["run"]        

        if check_if_failure:
            # saves result to outputs/re_run/...
            N_cont_doses = data["N_cont_doses"]
            DS_lim = data["DS_lim"]
            PP.re_run_cont(NDoses=NDoses, N_cont_doses=N_cont_doses, DS_lim=DS_lim, run_indices=[run])
        else:   
            grid_out = PP.re_run_grid(NDoses=NDoses, run_indices=[run], plot=plot)
            filename = f"./param_scan/outputs/failed/grid_output_run={run}.pickle"
            object_dump(filename, grid_out)






if __name__=="__main__":
    
    # N_dose_def = 11
    # N_cont_ds = 5
    
    N_dose_def = 101
    N_cont_ds = 300

    run_attrs = [
        # dict(run = 10,  DS_lim=[0.5,1.1], NDoses=N_dose_def, N_cont_doses=N_cont_ds), # Y
        # dict(run = 116, DS_lim=[0.2,1.0], NDoses=N_dose_def, N_cont_doses=N_cont_ds), # Y
        dict(run = 91,  DS_lim=[0.65,0.7], NDoses=N_dose_def, N_cont_doses=N_cont_ds), # N
        dict(run = 183, DS_lim=[0.3,0.6], NDoses=N_dose_def, N_cont_doses=N_cont_ds), # N
        dict(run = 227, DS_lim=[0.6,0.7], NDoses=N_dose_def, N_cont_doses=N_cont_ds), # N
        dict(run = 241, DS_lim=[0.3,0.6], NDoses=N_dose_def, N_cont_doses=N_cont_ds), # N
        dict(run = 295, DS_lim=[0.6,1.0], NDoses=N_dose_def, N_cont_doses=N_cont_ds), # N
        dict(run = 478, DS_lim=[0.2,0.6], NDoses=N_dose_def, N_cont_doses=N_cont_ds), # N

        # ESFY outperformed ERFB
        dict(run = 265, DS_lim=[0.4,0.8], NDoses=N_dose_def, N_cont_doses=N_cont_ds), # N
        ]
    
    
    
    
    if len(sys.argv)!=2:
        raise Exception("Supply one argument: the run index")
    index = int(sys.argv[1])
    run_attrs_use = [run_attrs[index]]



    # to check if failed
    main(run_attrs_use, plot=False, check_if_failure=True, print_it=True)
    
    # to re run failures via cluster
    # main(run_attrs_use, plot=False, check_if_failure=False, print_it=False)
    


