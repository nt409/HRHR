"""
Analyse the parameter scan output
"""
import pandas as pd

from .config import config
from .functions import PostProcess, get_par_str


def main():
    
    par_str = get_par_str(config)

    df = pd.read_csv(f"../outputs/csvs/PS_combined_output_{par_str}.csv")

    PP = PostProcess(df, par_str)

    PP.process_param_scan_df()

    PP.process_best_doses()

    PP.show_failed_runs()

    PP.check_max_EL_by_MS()

    PP.check_monotone_RFB()


if __name__=="__main__":
    main()