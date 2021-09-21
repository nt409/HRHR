"""
Re-run failed cases - try to find out why
"""

from param_scan.fns.config import config_rand
from param_scan.fns.post_process import PostProcess



    
PP = PostProcess(config_rand['par_str'])

df = PP.processed_df

print("FAILURES:")
print(df.loc[df["c_R_maxCont%"]<100, ["max_grid_EL", "c_R_maxContEL", "I_R_best_value"]])

print("ESFY 'outpeforms' ERFB:")
print(df.loc[df["c_R_maxContEL"]<df["c_E_maxContEL"], ["max_grid_EL", "c_E_maxContEL", "c_R_maxContEL", "I_R_best_value"]])

# * default 
# DS_lim = [0,2]

# run_indices=[117, 278, 290, 370]
run_indices=[17]

# run 290 was also outperformed by ESFY but fixed with denser grid
# also run 17 is outperformed by ESFY... want to get 12

run_IVT = False
run_cont = True

# 117 looks like a problem, others are ok?


if 278 in run_indices:
    if run_IVT:
        NDoses = 101
        PP.re_run_IVT(NDoses=NDoses, run_indices=[278])

    if run_cont:            
        NDoses = 101
        N_cont_doses = 101
        DS_lim = [0.30, 0.33]
        PP.re_run_cont(NDoses=NDoses, N_cont_doses=N_cont_doses, DS_lim=DS_lim, run_indices=[278])

if 290 in run_indices:
    if run_IVT:
        NDoses = 101
        PP.re_run_IVT(NDoses=NDoses, run_indices=[290])

    if run_cont:            
        NDoses = 101
        N_cont_doses = 101
        DS_lim = [0.1, 0.35]
        PP.re_run_cont(NDoses=NDoses, N_cont_doses=N_cont_doses, DS_lim=DS_lim, run_indices=[290])
        # PP.re_run_grid(NDoses=NDoses, run_indices=[290])

if 370 in run_indices:
    if run_IVT:
        NDoses = 101
        PP.re_run_IVT(NDoses=NDoses, run_indices=[370])
    
    if run_cont:            
        NDoses = 101
        N_cont_doses = 101
        DS_lim = [0.45, 0.65]
        PP.re_run_cont(NDoses=NDoses, N_cont_doses=N_cont_doses, DS_lim=DS_lim, run_indices=[370])
        # PP.re_run_grid(NDoses=NDoses, run_indices=[370])



if 117 in run_indices:
    if run_IVT:
        NDoses = 201
        PP.re_run_IVT(NDoses=NDoses, run_indices=[117])

    if run_cont:
        NDoses = 101
        N_cont_doses = 201
        DS_lim = [0.45, 0.56]
        PP.re_run_cont(NDoses=NDoses, N_cont_doses=N_cont_doses, DS_lim=DS_lim, run_indices=[117])


if 17 in run_indices:
    if run_IVT:
        NDoses = 201
        PP.re_run_IVT(NDoses=NDoses, run_indices=[17])

    if run_cont:
        NDoses = 101
        N_cont_doses = 201
        DS_lim = [0, 2]
        PP.re_run_cont(NDoses=NDoses, N_cont_doses=N_cont_doses, DS_lim=DS_lim, run_indices=[17])


