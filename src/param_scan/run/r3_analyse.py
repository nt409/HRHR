"""
Analyse the parameter scan output
"""

from ..fns.config import config_rand
from ..fns.post_process import PostProcess


def main(config):
    
    PP = PostProcess(config['folder_save'], config['par_str'])

    PP.get_maximum_along_contour_df()
    
    PP.analyse_max_contour_df()

    PP.analyse_failed()
    
    # PP.which_runs_worked_max_cont()

    # PP.re_run(NDoses=6, run_indices=list(range(1)))
    # PP.re_run(NDoses=41, run_indices=[706, 538, 548, 403, 490, 474, 367, 34])
    PP.re_run(NDoses=101, run_indices=[538])

    PP.check_high_or_low_dose()



if __name__=="__main__":
    main(config_rand)