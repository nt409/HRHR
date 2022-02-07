"""Now we use notebooks/006_dose_space_fig.ipynb to generate this figure, using
the functions below"""

import os
import pickle


from model.utils import object_dump
from model.strategy_arrays import EqualResFreqBreakdownArray, EqualSelectionArray
from model.simulator import RunGrid, RunSingleTactic
from model.config_classes import GridConfig, SingleConfig

from plotting.paper_figs import DoseSpaceOverview


class PretendPars:
    def __init__(self, primary_inoc) -> None:
        self.primary_inoc = primary_inoc
        self.fung_parms = None

    def get_single_conf(self, dose1, dose2):

        conf = SingleConfig(30, 10**(-7), 10**(-3),
                            dose1, dose1, dose2, dose2,
                            primary_inoculum=None)

        config_out = self._process_conf(conf)

        self.sing_conf = config_out

        return config_out

    def _process_conf(self, conf):

        conf.primary_inoculum = self.primary_inoc

        conf.bs_sex_prop = 0

        conf.load_saved = False

        conf.save = False

        conf.add_string()

        return conf


def get_contour_data(
    pars,
    output,
    n_contours,
    primary_inoc,
    load_saved,
):

    filename = (
        f"../outputs/saved_runs/contours_{n_contours}"
        f"_PI={primary_inoc['SR']:.8f}_{primary_inoc['RS']:.8f}"
        f"_{primary_inoc['RR']:.8f}.pickle"
    )

    if load_saved and os.path.isfile(filename):
        with open(filename, 'rb') as f:
            loaded = pickle.load(f)
        return loaded

    erfb_data = _get_this_contour_data(
        pars,
        EqualResFreqBreakdownArray(output),
        n_contours,
        primary_inoc,

    )

    esfy_data = _get_this_contour_data(
        pars,
        EqualSelectionArray(output),
        n_contours,
        primary_inoc,
    )

    run = dict(ERFB=erfb_data, ESFY=esfy_data)

    object_dump(filename, run)

    return run


def _get_this_contour_data(
    pars,
    contour_array_obj,
    n_contours,
    primary_inoc,
):

    cntrs = contour_array_obj.find_contours(pars, n_contours)

    FYs = []
    DSs = []
    xs = []
    ys = []

    for x, y in zip(cntrs['x'], cntrs['y']):
        conf_sing = SingleConfig(30, None, None, x, x, y, y)
        conf_sing.load_saved = False
        conf_sing.primary_inoculum = primary_inoc
        conf_sing.add_string()

        FY = RunSingleTactic().run(conf_sing).failure_year

        FYs.append(FY)
        DSs.append(x+y)
        xs.append(x)
        ys.append(y)

    out = dict(FY=FYs, DS=DSs, x=xs, y=ys)
    return out
