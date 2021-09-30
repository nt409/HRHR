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
    
    SREffect(df, n_its, filename)