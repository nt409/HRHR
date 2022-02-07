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


def combine(n_doses, n_its):
    df = pd.DataFrame()

    for alt_strat, index in itertools.product(["alt_12", "alt_21"], list(range(n_its))):

        filename = f"./alternation_scan/outputs/single/out_rand_{n_doses}_{alt_strat}_{index}.csv"

        new_df = pd.read_csv(filename)
        new_df["alt_strat"] = alt_strat

        df = pd.concat([df, new_df])

    # print(df.head(-10))
    # print(df.describe())

    print("WORKED RUNS; MIX AT LEAST AS GOOD")
    print(df[df["asex_alt"] <= df["asex_mix"]])
    print(df[df["sex_alt"] <= df["sex_mix"]])

    print("FAILED RUNS; ALT OUTPERFORMS MIX")
    print(df[df["asex_alt"] > df["asex_mix"]])
    print(df[df["sex_alt"] > df["sex_mix"]])

    print("NB 100 iterations but over alt_21 and alt_12")

    filename = f"./alternation_scan/outputs/combined/out_{n_its}_{n_doses}_{alt_strat}.csv"
    print(f"saving df to: {filename}")
    df.to_csv(filename, index=False)
