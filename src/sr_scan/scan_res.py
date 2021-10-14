import numpy as np
import pandas as pd
import sys

from model.simulator import RunGrid
from model.config_classes import GridConfig
from sr_scan.utils import get_sr_scan_params_res
from sr_scan.configs import config_res


def get_sr_scan_df_res(n_its, n_sex_props, n_doses, double_freq_factors, index):

    df = pd.DataFrame()
    
    sex_props = np.linspace(0, 1, n_sex_props)



    for dd in double_freq_factors:
        
        dfp = get_sr_scan_params_res(dd, index)



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


    dff_str = ",".join([str(ee) for ee in double_freq_factors])
    filename = f"./sr_scan/outputs/single/df_res_{n_its}_{n_sex_props}_{n_doses}_{dff_str}_{index}.csv"
    print(f"Saving df to: {filename}")
    df.to_csv(filename, index=False)
    
    return df











if __name__=="__main__":
    
    if len(sys.argv)!=2:
        raise Exception("Supply one argument: a run index")

    ind_dict = dict(index=int(sys.argv[1]))

    config_res = {**config_res, **ind_dict}
    
    df = get_sr_scan_df_res(**config_res)
        



