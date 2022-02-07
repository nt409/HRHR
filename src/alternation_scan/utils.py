"""Utility functions for alternation vs mixture scan"""

import numpy as np
import itertools
import pandas as pd

from model.config_classes import SingleConfig
from model.params import PARAMS
from model.simulator import RunSingleTactic


def random_alt_scan_params(index):
    """Randomly select rs, sr, rr and fcide params"""

    np.random.seed(index)
    is_invalid = True
    while is_invalid:
        rf1 = 10**(np.random.uniform(-10, -4))
        rf2 = 10**(np.random.uniform(-10, -4))
        rfd = 10**(np.random.uniform(-15, -6))

        om1 = np.random.uniform(0.4, 1)
        om2 = np.random.uniform(0.4, 1)
        thet1 = np.random.uniform(4, 12)
        thet2 = np.random.uniform(4, 12)
        delta_factor_1 = np.random.uniform(1/3, 3)
        delta_factor_2 = np.random.uniform(1/3, 3)

        primary_inoc = dict(RR=rfd, RS=rf1, SR=rf2, SS=1-rf1-rf2-rfd)

        fcide_parms = dict(
            omega_1=om1,
            omega_2=om2,
            theta_1=thet1,
            theta_2=thet2,
            delta_1=PARAMS.delta_1*delta_factor_1,
            delta_2=PARAMS.delta_2*delta_factor_2,
        )

        # check full dose fine:
        conf_a = SingleConfig(1, None, None, 1, 1, 1, 1)
        conf_a.load_saved = False
        conf_a.primary_inoculum = primary_inoc
        conf_a.add_string()

        yld = RunSingleTactic(fcide_parms).run(conf_a).yield_vec[0]

        if yld > 95:
            is_invalid = False

    return primary_inoc, fcide_parms


def combine_alt_scan_csvs(n_doses, n_its):
    df = pd.DataFrame()

    for index in range(n_its):

        filename = f"./alternation_scan/outputs/single/scan_rand_output_{n_doses}_{index}.csv"

        new_df = pd.read_csv(filename)

        df = pd.concat([df, new_df])

    filename = f"./alternation_scan/outputs/combined/out_{n_its}_{n_doses}.csv"
    print(f"saving df to: {filename}")
    df.to_csv(filename, index=False)
