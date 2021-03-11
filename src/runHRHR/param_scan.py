"""
Scan over: 
- IRFs
- fung asymptote
- fung decay rate
- partial res
- levels of SR
"""
import numpy as np

from utils.functions import process_param_scan_df, process_best_doses, \
    show_failed_runs, run_param_scan, check_max_EL_by_MS, check_monotone_RFB


to_scan_over = dict(
    RFS1 = 10**(-5),
    RFS2 = [10**(k) for k in np.linspace(-8,-1, 2)],
    RFD = [10**(k) for k in np.linspace(-15,-1, 2)],

    asym1 = np.linspace(0.5, 1, 4),
    asym2 = np.linspace(0.5, 1, 4),

    dec_rate1 = [(1.11*10**(-2))*(2**j) for j in [-2,0,2]],
    dec_rate2 = [(1.11*10**(-2))*(2**j) for j in [-2,-1,0,1,2]],

    SR = list(np.linspace(0,1,5)) + [0.9, 0.95, 0.99],
    )

df = run_param_scan(to_scan_over)

best_doses = process_param_scan_df(df)

delt_worked, all_runs, runs_worked, runs_failed = process_best_doses(best_doses)

show_failed_runs(df, runs_failed)

max_EL_by_MS = check_max_EL_by_MS(df)

mon_RFB = check_monotone_RFB(df)