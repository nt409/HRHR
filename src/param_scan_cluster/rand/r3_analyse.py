"""
Analyse the parameter scan output
"""
import pandas as pd

from ..config import config_rand
from ..functions import PostProcess, get_PS_rand_str


def main(config):
    
    par_str = get_PS_rand_str(config)

    df = pd.read_csv(f"param_scan_cluster/outputs/rand/combined/output_{par_str}.csv")

    PP = PostProcess(df, par_str)

    if config['contour_type']=="FYY":
        # for i in range(10):
        #     PP.testing_RFB(i)

        # PP.process_best_doses()

        PP.check_max_EL_by_contour("rand")

        PP.check_monotone_RFB("rand")
        
        PP.which_runs_worked_monotone_RFB(print_=True)

        # PP.plot_df(run=0)
        
        # PP.maxEL_by_contour()
        
        NDOSES = 15

        PP.re_run_failures(NDOSES, failed_run_indices=list(range(6)))
    
    else:
        PP.get_maximum_along_contour_df()
        
        PP.analyse_max_contour_df()

        PP.analyse_failed()
        
        PP.which_runs_worked_max_cont()

        PP.re_run_failures(NDoses=11, failed_run_indices=None)
        
        # PP.check_high_or_low_dose()

if __name__=="__main__":
    main(config_rand)