"""Alternation vs mixtures scan

Called by alternation_scan.scan.submit

- Randomly select rs, sr, rr and fcide params
- Check best mixture vs best alternation 12 & 21

"""

import copy
import sys
import numpy as np
import pandas as pd


from model.simulator import RunGrid
from model.config_classes import GridConfig
from alternation_scan.utils import random_alt_scan_params


def main(n_doses, index, load_saved):

    primary_inoc, fcide_parms = random_alt_scan_params(index)

    model_run = RunGrid(fcide_parms)

    # MIXTURE
    # asex, mixture
    cnfg_asex_mix = GridConfig(40, None, None, n_doses)
    cnfg_asex_mix.load_saved = load_saved

    cnfg_asex_mix.primary_inoculum = primary_inoc
    cnfg_asex_mix.add_string()

    asex_mix = model_run.run(cnfg_asex_mix)

    # sex, mixture
    cnfg_sex_mix = copy.copy(cnfg_asex_mix)
    cnfg_sex_mix.bs_sex_prop = 1
    cnfg_sex_mix.add_string()

    sex_mix = model_run.run(cnfg_sex_mix)

    # ALTERNATION
    # asex, alt 21
    cnfg_asex_alt_21 = copy.copy(cnfg_asex_mix)
    cnfg_asex_alt_21.strategy = "alt_21"
    cnfg_asex_alt_21.add_string()

    asex_alt_21 = model_run.run(cnfg_asex_alt_21)

    # sex, alt 21
    cnfg_sex_alt_21 = copy.copy(cnfg_asex_alt_21)
    cnfg_sex_alt_21.bs_sex_prop = 1
    cnfg_sex_alt_21.add_string()

    sex_alt_21 = model_run.run(cnfg_sex_alt_21)

    # ALTERNATION
    # asex, alt 12
    cnfg_asex_alt_12 = copy.copy(cnfg_asex_mix)
    cnfg_asex_alt_12.strategy = "alt_12"
    cnfg_asex_alt_12.add_string()

    asex_alt_12 = model_run.run(cnfg_asex_alt_12)

    # sex, alt 12
    cnfg_sex_alt_12 = copy.copy(cnfg_asex_alt_12)
    cnfg_sex_alt_12.bs_sex_prop = 1
    cnfg_sex_alt_12.add_string()

    sex_alt_12 = model_run.run(cnfg_sex_alt_12)

    # ADD TO DATAFRAME
    data = dict(
        RR=primary_inoc["RR"],
        RS=primary_inoc["RS"],
        SR=primary_inoc["SR"],

        omega1=fcide_parms["omega_1"],
        omega2=fcide_parms["omega_2"],
        theta1=fcide_parms["theta_1"],
        theta2=fcide_parms["theta_2"],
        delta1=fcide_parms["delta_1"],
        delta2=fcide_parms["delta_2"],

        asex_alt_21=np.amax(asex_alt_21.FY),
        sex_alt_21=np.amax(sex_alt_21.FY),

        asex_alt_12=np.amax(asex_alt_12.FY),
        sex_alt_12=np.amax(sex_alt_12.FY),

        asex_mix=np.amax(asex_mix.FY),
        sex_mix=np.amax(sex_mix.FY),
    )

    df = pd.DataFrame([data])

    filename = f"./alternation_scan/outputs/single/scan_rand_output_{n_doses}_{index}.csv"
    print(f"saving df to: {filename}")
    df.to_csv(filename)

    return df


if __name__ == "__main__":

    if len(sys.argv) != 2:
        raise Exception("Supply one argument: the run index")

    index = int(sys.argv[1])

    df = main(n_doses=51, index=index, load_saved=False)
