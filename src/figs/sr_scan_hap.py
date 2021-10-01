from sr_hap.analyse import get_haploid_outputs
import pandas as pd
import pandas as pd

from plotting.paper_figs import SREffect



if __name__=="__main__":
    
    n_variants = 3
    n_sex_props = 11
    n_doses = 21

    n_trcs_per_fig = 9
    
    # n_variants = 5
    # n_sex_props = 2
    # n_doses = 5
    

    double_freq_factor_lowest = 1e-4
    dff_str = ",".join([str(ee) for ee in double_freq_factor_lowest])
    filename = f"./sr_hap/outputs/combined/df_{n_variants}_{n_sex_props}_{n_doses}_{dff_str}.csv"
    df = pd.read_csv(filename)    

    indices = [6, 10, 24]


    outputs_l = get_haploid_outputs(n_variants, n_doses, double_freq_factor_lowest, indices[0], [0,0.5,1])
    outputs_m = get_haploid_outputs(n_variants, n_doses, double_freq_factor_lowest, indices[1], [0,1])
    outputs_h = get_haploid_outputs(n_variants, n_doses, double_freq_factor_lowest, indices[2], [0,1])


    filename = f"../outputs/figures/paper_figs/sr_effect_{n_variants}_{n_sex_props}_{n_doses}_{dff_str}.png"
    SREffect(df, outputs_l, outputs_m, outputs_h, n_trcs_per_fig, indices, filename)