import numpy as np
import pandas as pd
import sys

from model.simulator import RunGrid, RunSingleTactic
from model.config_classes import GridConfig
from sr_scan.utils import get_sr_scan_params_app
from sr_scan.configs import config_res


def get_sr_scan_df_res(n_its, n_sex_props, n_doses, double_freq_factors, index):

    
    sex_props = np.linspace(0, 1, n_sex_props)

    df = pd.DataFrame()

    for dd in double_freq_factors:
        tmp = pd.DataFrame()

        dfp, fcide_pars, conf = get_sr_scan_params_app(dd, index)
        
        conf.load_saved = False
        



        data = dfp.iloc[int(0),:]


        conf.load_saved = False


        for bs in sex_props:

            conf.bs_sex_prop = bs
            conf.add_string()

            out = RunSingleTactic(fcide_pars).run(conf)


            data.loc["maxDoseEL"] = out.failure_year
            data.loc["bs_sex_prop"] = bs
            data.loc["run"] = index
            
            tmp = tmp.append(data, ignore_index=True)

        increasing = np.float(tmp.loc[((tmp["bs_sex_prop"]==1) & (tmp["double_freq_factor"]==dd)), "maxDoseEL"])> np.float(tmp.loc[((tmp["bs_sex_prop"]==0) & (tmp["double_freq_factor"]==dd)), "maxDoseEL"])
        const      = np.float(tmp.loc[((tmp["bs_sex_prop"]==1) & (tmp["double_freq_factor"]==dd)), "maxDoseEL"])==np.float(tmp.loc[((tmp["bs_sex_prop"]==0) & (tmp["double_freq_factor"]==dd)), "maxDoseEL"])
        decreasing = np.float(tmp.loc[((tmp["bs_sex_prop"]==1) & (tmp["double_freq_factor"]==dd)), "maxDoseEL"])< np.float(tmp.loc[((tmp["bs_sex_prop"]==0) & (tmp["double_freq_factor"]==dd)), "maxDoseEL"])

        this_bit_worked = ((increasing and dd>1) or (const and dd==1) or (decreasing and dd<1))
        worked = True if this_bit_worked else False

        tmp["worked"] = worked

        if any(tmp["worked"]):
            df = pd.concat([df, tmp])
    
    if any(df["worked"]):
        dff_str = ",".join([str(ee) for ee in double_freq_factors])
        filename = f"./sr_scan/outputs/single/df_app_{n_its}_{n_sex_props}_{n_doses}_{dff_str}_{index}.csv"
        print(f"Saving df to: {filename}")
        df.to_csv(filename, index=False)
    
    return df











if __name__=="__main__":
    
    df = pd.DataFrame()

    for ind in range(10):

        ind_dict = dict(index=ind)

        config_res = {**config_res, **ind_dict}
        
        tmp = get_sr_scan_df_res(**config_res)

        if any(tmp["worked"]):
            df = pd.concat([df, tmp])
    
    filename = f"./sr_scan/outputs/combined/df_app.csv"
    print(f"Saving df to: {filename}")
    df.to_csv(filename, index=False)
        



