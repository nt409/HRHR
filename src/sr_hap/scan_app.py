
import numpy as np
import pandas as pd
import sys

from model.simulator import RunGrid
from model.config_classes import GridConfig
from .utils import get_sr_scan_params_app



def get_sr_scan_df_app(n_variants, n_sex_props, n_doses, double_freq_factor_lowest, index):
    
    dfp = get_sr_scan_params_app(n_variants, double_freq_factor_lowest, index)

    sex_props = np.linspace(0, 1, n_sex_props)

    df = pd.DataFrame()



    data = dfp.iloc[int(0),:]

    conf = GridConfig(30, None, None, n_doses)
    
    conf.primary_inoculum = dict(
        RR = data["RR"],
        RS = data["RS"],
        SR = data["SR"],
        SS = data["SS"],
        )

    fcide_pars = dict(
        omega_1 = data["omega_1"],
        omega_2 = data["omega_2"],
        theta_1 = data["theta_1"],
        theta_2 = data["theta_2"],
        delta_1 = data["delta_1"],
        delta_2 = data["delta_2"],
        )

    conf.load_saved = False

    for bs in sex_props:

        conf.bs_sex_prop = bs
        conf.add_string()

        output = RunGrid(fcide_pars).run(conf)

        data["maxEL"] = np.amax(output.FY)
        data["bs_sex_prop"] = bs
        data["run"] = index
        
        df = df.append(data, ignore_index=True)



    filename = f"./sr_hap/outputs/single/df_app_{n_variants}_{n_sex_props}_{n_doses}_{double_freq_factor_lowest}_{index}.csv"
    print(f"Saving df to: {filename}")
    df.to_csv(filename, index=False)
    
    return df








if __name__=="__main__":
    n_variants = 3
    n_sex_ps = 11
    n_doses = 21
    
    if len(sys.argv)!=2:
        raise Exception("Supply one argument: a run index")

    index = int(sys.argv[1])
    
    double_freq_factor_lowest = 1e-4
    df = get_sr_scan_df_app(n_variants, n_sex_ps, n_doses, double_freq_factor_lowest, index)
        



