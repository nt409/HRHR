import itertools
import numpy as np
import pandas as pd
import os
from tqdm import tqdm

from model.simulator import RunSingleTactic
from model.config_classes import SingleConfig


def get_sf_ratio_data(n_doses, fcide_pars=None, save=True):
    fp = fcide_pars
    fcide_str = "default" if fp is None else f"{fp['theta_1']}_{fp['omega_1']}"
    filename = f"./sr_scan/outputs/selection_ratio_{n_doses}_{fcide_str}.csv"
    
    if os.path.isfile(filename):
        loaded_df = pd.read_csv(filename)
        return loaded_df



    doses = np.linspace(0, 1, n_doses)
    ifs = [1e-5, 1e-1]

    df = pd.DataFrame()

    for rf, d in tqdm(itertools.product(ifs, doses)):
        config_sing = SingleConfig(1, None, None, d, d, d, d)
        config_sing.primary_inoculum = dict(RR=rf**2, RS=rf, SR=rf, SS=1-2*rf-rf**2)
        config_sing.load_saved = False

        config_sing.ws_sex_prop = 0
        config_sing.bs_sex_prop = 1
        config_sing.add_string()

        output = RunSingleTactic(fcide_pars).run(config_sing)

        sf_ratio_sing = output.start_freqs['SR'][1]/output.start_freqs['SR'][0]
        sf_ratio_doub = output.start_freqs['RR'][1]/output.start_freqs['RR'][0]

        data = dict(d=d, 
                    sf_ratio_sing=sf_ratio_sing,
                    sf_ratio_doub=sf_ratio_doub,
                    SR=output.start_freqs['SR'][0],
                    RR=output.start_freqs['RR'][0],
                    yld=output.yield_vec[0])

        df = df.append(data, ignore_index=True)

    df = df.loc[df["yld"]>95]

    if save:
        print(f"Saving df to: {filename}")
        df.to_csv(filename)
    
    return df



if __name__=="__main__":
    get_sf_ratio_data(5)
