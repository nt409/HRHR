"""
Analyse the parameter scan output
"""
import pandas as pd

from .config import config
from .functions import PostProcess, get_PS_grid_str


def main(config):
    
    par_str = get_PS_grid_str(config)

    df = pd.read_csv(f"param_scan_cluster/outputs/PS_combined_grid_output_{par_str}.csv")

    PP = PostProcess(df, par_str)

    PP.process_param_scan_df()

    PP.process_best_doses()

    PP.show_failed_runs()

    PP.check_max_EL_by_MS("grid")

    PP.check_monotone_RFB("grid")


if __name__=="__main__":
    main(config)