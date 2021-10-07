import pandas as pd

from sr_hap.configs import config_res
from plotting.paper_figs import SREffectResults, SREffectResults2



if __name__=="__main__":

    n_its = config_res["n_its"]
    n_sex_props = config_res["n_sex_props"]
    n_doses = config_res["n_doses"]
    double_freq_factors = config_res["double_freq_factors"]
    

    dff_str = ",".join([str(ee) for ee in double_freq_factors])
    filename = f"./sr_hap/outputs/combined/df_res_{n_its}_{n_sex_props}_{n_doses}_{dff_str}.csv"
    df = pd.read_csv(filename)


    filename = f"../outputs/figures/paper_figs/sr_effect_res_{n_its}_{n_sex_props}_{n_doses}_{dff_str}.png"
    SREffectResults2(df, double_freq_factors, filename)
