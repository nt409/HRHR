"""
Analyse the parameter scan output, and get table of results for paper
"""

from param_scan.fns.config import config_rand
from param_scan.fns.post_process import PostProcess


def main(config):
    
    PP = PostProcess(config['par_str'])

    mcdf = PP.processed_df
    
    print(mcdf.loc[:, ["ERFB_diff_from_opt"]].value_counts())
    print(mcdf.loc[:, ["ESFY_diff_from_opt"]].value_counts())
    
    results_df = PP.analyse_processed_df()
    print(results_df)




if __name__=="__main__":
    main(config_rand)