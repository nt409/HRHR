"""
Analyse the parameter scan output, and get table of results for paper
"""

from param_scan.fns.config import config_rand
from param_scan.fns.post_process import PostProcess


def main(config):

    PP = PostProcess(config['par_str'])

    dfp = PP.processed_df

    # print(dfp.loc[dfp["c_R_maxContEL"] -4 >= dfp["c_E_maxContEL"], ["c_R_maxContEL", "c_E_maxContEL"]])

    print("\nMean max EL if: (i) general; (ii) failed by at least 4:")
    print(dfp["c_R_maxContEL"].mean())
    print(dfp.loc[dfp["c_R_maxContEL"] - 4 >= dfp["c_E_maxContEL"],
          ["c_R_maxContEL", "c_E_maxContEL"]].mean())

    print("\nThese runs failed:")
    print(dfp.loc[dfp["c_R_maxContEL"] < dfp["max_grid_EL"], ["run"]])

    print("\nThese runs had ESFY outperformed ERFB:")
    print(dfp.loc[dfp["c_R_maxContEL"] < dfp["c_E_maxContEL"], [
          "run", "c_R_maxContEL", "c_E_maxContEL", "max_grid_EL"]])

    print("\nThis is how different the outcomes were:")
    print(dfp.loc[:, ["ERFB_diff_from_opt"]].value_counts())
    print(dfp.loc[:, ["ESFY_diff_from_opt"]].value_counts())
    print(dfp.loc[:, ["ESFYL_diff_from_opt"]].value_counts())

    print("\nESFY/ESFYL wrong by at least two years:")
    print(100*sum(dfp.loc[dfp["ESFY_diff_from_opt"]
          <= -2, ["ESFY_diff_from_opt"]].count())/500)
    print(100*sum(dfp.loc[dfp["ESFYL_diff_from_opt"]
          <= -2, ["ESFYL_diff_from_opt"]].count())/500)

    results_df = PP.analyse_processed_df()
    print("\nThis is the summary table for the paper:")
    print(results_df)


if __name__ == "__main__":
    main(config_rand)
