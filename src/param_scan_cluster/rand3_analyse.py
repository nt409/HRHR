"""
Analyse the parameter scan output
"""
import pandas as pd

from .config import config_rand
from .functions import PostProcess, get_PS_rand_str


def main(config):
    
    par_str = get_PS_rand_str(config)

    df = pd.read_csv(f"param_scan_cluster/outputs/PS_combined_rand_output_{par_str}.csv")
    # print(df.run.unique())
    # exit()

    PP = PostProcess(df, par_str)

    PP.process_param_scan_df()

    PP.process_best_doses()

    PP.show_failed_runs()

    PP.check_max_EL_by_MS("rand")

    PP.check_monotone_RFB("rand")


if __name__=="__main__":
    main(config_rand)