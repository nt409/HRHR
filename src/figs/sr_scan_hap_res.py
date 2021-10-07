import pandas as pd

from sr_hap.utils import get_haploid_outputs_res
from plotting.paper_figs import SREffectResults



if __name__=="__main__":

    n_its = 10
    n_sex_props = 11
    n_doses = 21
    

    double_freq_factors = [1e-5, 1, 1e5]
    dff_str = ",".join([str(ee) for ee in double_freq_factors])
    filename = f"./sr_hap/outputs/combined/df_res_{n_its}_{n_sex_props}_{n_doses}_{dff_str}.csv"
    df = pd.read_csv(filename)

    indices = [4, 13, 22]

    indices_use = list(range(27))

    outputs_l = get_haploid_outputs_res(n_its, n_doses, double_freq_factors, indices[0], [0,0.2,1])
    outputs_m = get_haploid_outputs_res(n_its, n_doses, double_freq_factors, indices[1], [0,1])
    outputs_h = get_haploid_outputs_res(n_its, n_doses, double_freq_factors, indices[2], [0,1])


    filename = f"../outputs/figures/paper_figs/sr_effect_res_{n_its}_{n_sex_props}_{n_doses}_{dff_str}.png"

    SREffectResults(df, outputs_l, outputs_m, outputs_h, n_its, filename)
