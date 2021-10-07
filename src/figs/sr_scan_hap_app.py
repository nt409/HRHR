import pandas as pd

from sr_hap.utils import get_haploid_outputs_app
from plotting.paper_figs import SREffectAppendix



if __name__=="__main__":

    n_variants = 3
    n_sex_props = 11
    n_doses = 21
    
    n_trcs_per_fig = 9
    

    double_freq_factor_lowest = 1e-4
    filename = f"./sr_hap/outputs/combined/df_app_{n_variants}_{n_sex_props}_{n_doses}_{double_freq_factor_lowest}.csv"
    df = pd.read_csv(filename)

    indices = [4, 13, 22]

    indices_use = list(range(27))

    outputs_l = get_haploid_outputs_app(n_variants, n_doses, double_freq_factor_lowest, indices[0], [0,0.2,1])
    outputs_m = get_haploid_outputs_app(n_variants, n_doses, double_freq_factor_lowest, indices[1], [0,1])
    outputs_h = get_haploid_outputs_app(n_variants, n_doses, double_freq_factor_lowest, indices[2], [0,1])


    filename = f"../outputs/figures/paper_figs/sr_effect_app_{n_variants}_{n_sex_props}_{n_doses}_{double_freq_factor_lowest}.png"
    SREffectAppendix(df, outputs_l, outputs_m, outputs_h, n_trcs_per_fig, indices, indices_use, filename)