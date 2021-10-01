from sr_hap.analyse import get_haploid_outputs
import pandas as pd
import pandas as pd

from plotting.paper_figs import SREffect



if __name__=="__main__":
    n_its = 10
    n_sex_props = 11
    n_doses = 21
    
    # n_its = 5
    # n_sex_props = 2
    # n_doses = 5
    

    double_freq_factors = [1e-5, 1, 1e5]
    dff_str = ",".join([str(ee) for ee in double_freq_factors])
    filename = f"./sr_hap/outputs/combined/df_{n_its}_{n_sex_props}_{n_doses}_{dff_str}.csv"
    df = pd.read_csv(filename)    


    outputs_l = get_haploid_outputs(n_its, n_doses, double_freq_factors, 6, [0,0.5,1])
    outputs_m = get_haploid_outputs(n_its, n_doses, double_freq_factors, 10, [0,1])
    outputs_h = get_haploid_outputs(n_its, n_doses, double_freq_factors, 24, [0,1])
    

    filename = f"../outputs/figures/paper_figs/sr_effect_{n_its}_{n_sex_props}_{n_doses}_{dff_str}.png"
    SREffect(df, outputs_l, outputs_m, outputs_h, n_its, filename)