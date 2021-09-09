"""
Re-run failed cases - try to find out why
"""

from param_scan.fns.config import config_rand
from param_scan.fns.post_process import PostProcess


def main(config):
    
    PP = PostProcess(config['par_str'])

    df = PP.processed_df

    print("FAILURES:")
    print(df.loc[df["c_R_maxCont%"]<100, ["max_grid_EL", "c_R_maxContEL", "I_R_best_value"]])


    # NDoses = 51
    # N_cont_doses = 3
    # run_indices=[117]
    
    NDoses = 51
    N_cont_doses = 101
    run_indices=[117, 278, 290, 370]


    PP.re_run_cont(NDoses=NDoses, N_cont_doses=N_cont_doses, run_indices=run_indices)
    
    # PP.re_run_grid(NDoses=NDoses, run_indices=run_indices)





if __name__=="__main__":
    main(config_rand)