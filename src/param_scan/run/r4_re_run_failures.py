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



# run_IVT = False
# run_cont = True

N_dose_def = 51
N_cont_ds = 100

run_attrs = [
    dict(run = 80,  DS_lim=[0.25,0.38], NDoses=N_dose_def, N_cont_doses=N_cont_ds), # Y,L -- ESFY one too
    dict(run = 279, DS_lim=[0.2,0.28], NDoses=N_dose_def, N_cont_doses=N_cont_ds),  # Y,L -- ESFY one too
    dict(run = 46,  DS_lim=[0.34,0.4], NDoses=N_dose_def, N_cont_doses=N_cont_ds),  # Y,L
    dict(run = 135, DS_lim=[0.3,0.38], NDoses=N_dose_def, N_cont_doses=N_cont_ds),  # Y,L
    dict(run = 165, DS_lim=[0.3,0.36], NDoses=N_dose_def, N_cont_doses=N_cont_ds),  # Y,L
    dict(run = 471, DS_lim=[0.2,0.3], NDoses=N_dose_def, N_cont_doses=N_cont_ds),   # Y,L

    # dict(run = 496, DS_lim=[0.43,0.45], NDoses=N_dose_def, N_cont_doses=150),   # N
    # dict(run = 114, DS_lim=[0.28,0.36], NDoses=N_dose_def, N_cont_doses=150), # N
    # dict(run = 143, DS_lim=[0.21,0.25], NDoses=N_dose_def, N_cont_doses=150),   # N
    
    dict(run = 328, DS_lim=[1.8,1.9],   NDoses=N_dose_def, N_cont_doses=150),   # N - probs not
    ]


for data in run_attrs:
    NDoses = data["NDoses"]
    N_cont_doses = data["N_cont_doses"]
    DS_lim = data["DS_lim"]
    run = data["run"]
    
    # if run_IVT:
    #     NDoses = 201
    #     PP.re_run_IVT(NDoses=NDoses, run_indices=[17])

    PP.re_run_grid(NDoses=NDoses, run_indices=[run])
    # PP.re_run_cont(NDoses=NDoses, N_cont_doses=N_cont_doses, DS_lim=DS_lim, run_indices=[run])


