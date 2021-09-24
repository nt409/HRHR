import copy
import itertools
import numpy as np
import pandas as pd
from tqdm import tqdm

from model.params import PARAMS
from model.simulator import RunGrid
from model.config_classes import GridConfig


def main(n_doses, n_its, within_season, alt_strat):

    rfs = np.logspace(-10, -4, n_its)
    omegas = np.linspace(0.4, 1, n_its)
    thetas = np.linspace(4,  12, n_its)

    df = pd.DataFrame()

    for rf1, rf2, om1, om2, thet1, thet2 in tqdm(itertools.product(rfs, [1e-5], omegas, [1], thetas, [9.6])):

        primary_inoc = dict(RR=rf1*rf2, RS=rf1, SR=rf2, SS=1-rf1-rf2-rf1*rf2)

        fcide_parms = dict(
            omega_1 = om1,
            omega_2 = om2,
            theta_1 = thet1,
            theta_2 = thet2,
            delta_1 = PARAMS.delta_1,
            delta_2 = PARAMS.delta_2,
            )

        model_run = RunGrid(fcide_parms)
        
        # MIXTURE
        # asex, mixture
        conf_grid_a_mix = GridConfig(40, None, None, n_doses)
        conf_grid_a_mix.load_saved = False
        
        if within_season:
            conf_grid_a_mix.ws_sex_prop = 1
        # else do nothing - default is 0

        conf_grid_a_mix.primary_inoculum = primary_inoc
        conf_grid_a_mix.add_string()

        asex_mix = model_run.run(conf_grid_a_mix)

        # sex, mixture
        conf_grid_s_mix = copy.copy(conf_grid_a_mix)
        conf_grid_s_mix.bs_sex_prop = 1
        conf_grid_s_mix.add_string()

        sex_mix = model_run.run(conf_grid_s_mix)

        # ALTERNATION
        # asex, alt
        conf_grid_a_alt = copy.copy(conf_grid_a_mix)
        conf_grid_a_alt.strategy = alt_strat
        conf_grid_a_alt.add_string()

        asex_alt = model_run.run(conf_grid_a_alt)

        # sex, alt
        conf_grid_s_alt = copy.copy(conf_grid_a_alt)
        conf_grid_s_alt.bs_sex_prop = 1
        conf_grid_s_alt.add_string()
        
        sex_alt = model_run.run(conf_grid_s_alt)

        # ADD TO DATAFRAME

        data = dict(
            rf1 =  rf1,
            rf2 = rf2,
            omega1 =  om1,
            omega2 =  om2,
            theta1 =  thet1,
            theta2 = thet2,
            asex_alt = np.amax(asex_alt.FY),
            asex_mix = np.amax(asex_mix.FY),
            sex_alt = np.amax(sex_alt.FY),
            sex_mix = np.amax(sex_mix.FY),
            )
        
        df = df.append(data, ignore_index=True)

    filename = f"./alt_scan/output_{n_its}_{n_doses}_{within_season}_{alt_strat}.csv"
    print(f"saving df to: {filename}")
    df.to_csv(filename)

    return df




if __name__=="__main__":

    # df = main(n_doses=21, n_its=5, within_season=False, alt_strat="alt_21")
    # df = main(n_doses=21, n_its=5, within_season=True, alt_strat="alt_21")
    # df = main(n_doses=21, n_its=5, within_season=False, alt_strat="alt_12")
    df = main(n_doses=21, n_its=5, within_season=True, alt_strat="alt_12")


