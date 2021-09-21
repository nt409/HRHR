import numpy as np
import pandas as pd
import itertools
import os
from tqdm import tqdm

from model.simulator import RunSingleTactic
from model.config_classes import SingleConfig


def get_sr_grid_df(n_doses, n_sex_props, fcide_pars=None, save=True):
    fp = fcide_pars
    fcide_str = "default" if fp is None else f"{fp['theta_1']}_{fp['omega_1']}"
    filename = f"./sr_scan/outputs/sr_grid_{n_doses}_{n_sex_props}_{fcide_str}.csv"
    
    if os.path.isfile(filename):
    # if False:
        loaded_df = pd.read_csv(filename)
        return loaded_df



    doses = np.linspace(0, 1, n_doses)
    sex_props = np.linspace(0, 1, n_sex_props)

    df = pd.DataFrame()

    for d, ws, bs in tqdm(itertools.product(doses, sex_props, sex_props)):
        config_sing = SingleConfig(30, None, None, d, d, d, d)
        rf = 1e-5
        config_sing.primary_inoculum = dict(RR=rf**2, RS=rf, SR=rf, SS=1-2*rf-rf**2)
        config_sing.load_saved = False

        config_sing.ws_sex_prop = ws
        config_sing.bs_sex_prop = bs
        config_sing.add_string()

        output = RunSingleTactic(fcide_pars).run(config_sing)

        data = dict(ws=ws, bs=bs, d=d, EL=output.failure_year)

        df = df.append(data, ignore_index=True)
    
    res = process_sr_grid(df)
    
    if save:
        print(f"Saving df to: {filename}")
        res = res.reset_index()
        res.to_csv(filename)
    
    return res











def process_sr_grid(df):
    grouped = df.groupby(['bs', 'ws'])
    
    low_d = grouped.apply(get_low_dose_EL)

    full_d = grouped.apply(get_full_dose_EL)

    left = low_d.set_index(['bs', 'ws'])
    right = full_d.set_index(['bs', 'ws'])
    
    res = left.join(right, lsuffix="_low", rsuffix="_full")

    res['Z_metric'] = res['EL_full']/(res['EL_low'] + res['EL_full'])

    return res




def get_low_dose_EL(df):
    valid = df.loc[df['EL']>1]

    min_d = valid.d.min()
    limit = min_d + (1/3)*(1-min_d)

    out = valid.loc[valid.d<=limit]
    out = out.loc[out.EL==out.EL.max(), :].iloc[0]

    out['min_valid_d'] = min_d
    out['thresh_d'] = limit
    
    return out

def get_full_dose_EL(df):
    fd = df.loc[df.d==1, ['bs', 'ws', 'd', 'EL']]
    return fd
