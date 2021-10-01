from math import log10
from model.simulator import RunGrid
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px


from model.config_classes import GridConfig
from sr_hap.scan import get_sr_scan_params


def get_haploid_outputs(n_variants, n_doses, double_freq_factor_lowest, index, bss):
    outputs = []

    for bs in bss:

        dfp = get_sr_scan_params(n_variants, double_freq_factor_lowest, index)

        output = check_run(n_doses, bs, dfp)

        outputs.append(output)
    return outputs






def check_run(n_doses, bs, dfp):

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

    conf.load_saved = True

    conf.bs_sex_prop = bs
    conf.add_string()

    output = RunGrid(fcide_pars).run(conf)
    return output

