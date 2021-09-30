import pandas as pd

from plotting.paper_figs import SREffect



if __name__=="__main__":
    n_its = 10
    n_sex_ps = 6
    n_doses = 21
    
    n_its = 5
    n_sex_ps = 2
    n_doses = 5

    res = pd.DataFrame()

    double_freqs = [-20, -12, -6]
    
    for double_freq in double_freqs:
        
        n_its = 10
        n_sex_ps = 6
        n_doses = 21
        filename = f"./sr_hap/outputs/df_{n_its}_{n_sex_ps}_{n_doses}_{double_freq}.csv"

        df = pd.read_csv(filename)
        res = pd.concat([res, df])
    
    dfr_str = ",".join([str(ee) for ee in double_freqs])
    filename = f"../outputs/figures/paper_figs/sr_effect_{n_its}_{n_sex_ps}_{n_doses}_{dfr_str}.png"
    
    SREffect(res, double_freqs, filename)