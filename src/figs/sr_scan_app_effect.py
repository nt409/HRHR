import copy
import pandas as pd

from model.simulator import RunSingleTactic
from plotting.paper_figs import SREffectAppendix
from sr_scan.configs import config_app
from sr_scan.utils import get_sr_scan_params_app


def main(dd, index, sex_props):

    _, fcide_pars, conf = get_sr_scan_params_app(dd, index)
    conf.load_saved = False

    out = []

    for bs in sex_props:
        conf.bs_sex_prop = bs
        conf.add_string()
        data = RunSingleTactic(fcide_pars).run(conf)
        out.append(data)

    conf_dict = copy.copy(vars(conf))

    print(pd.DataFrame([fcide_pars]).iloc[0])
    print(conf_dict['primary_inoculum'])

    return out


if __name__ == "__main__":

    sex_props = [0, 0.5, 1]

    outputs_l = main(config_app["double_freq_factors"]
                     [0], index=3, sex_props=sex_props)
    outputs_m = main(config_app["double_freq_factors"]
                     [1], index=4, sex_props=sex_props)
    outputs_h = main(config_app["double_freq_factors"]
                     [2], index=9, sex_props=sex_props)

    filename = f"../outputs/figures/paper_figs/sr_effect_app.png"
    SREffectAppendix(outputs_l, outputs_m, outputs_h, sex_props, filename)
