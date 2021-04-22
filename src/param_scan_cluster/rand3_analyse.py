"""
Analyse the parameter scan output
"""
import pandas as pd

from .config import config_rand
from .functions import PostProcess, get_PS_rand_str


def main(config):
    
    par_str = get_PS_rand_str(config)

    df = pd.read_csv(f"param_scan_cluster/outputs/PS_combined_rand_output_{par_str}.csv")

    PP = PostProcess(df, par_str)

    # PP.process_best_doses()

    PP.check_max_EL_by_MS("rand")

    PP.check_monotone_RFB("rand")
    
    PP.which_runs_worked(print_=True)
    
    PP.re_run_failures(config["NDoses"], failed_run_indices=list(range(3,6)))


if __name__=="__main__":
    main(config_rand)