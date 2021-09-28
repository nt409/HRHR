from model.utils import object_dump
import os
import pickle
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






def get_contour_data(n_contours, rf1, rf2, load_saved):
    
    filename = f'../outputs/saved_runs/contours_{n_contours}_{rf1}_{rf2}.pickle'

    if load_saved and os.path.isfile(filename):
        with open(filename, 'rb') as f:
            loaded = pickle.load(f)
        return loaded
    
    d1 = get_this_contour_data(EqualResFreqBreakdownArray(output), n_contours, rf1, rf2)
    d2 = get_this_contour_data(EqualSelectionArray(output), n_contours, rf1, rf2)
    run = dict(ERFB=d1, ESFY=d2)
    object_dump(filename, run)

    return run




def get_this_contour_data(z, n_contours, rf1, rf2):
    cntrs = z.find_contours(pars, n_contours)

    FYs = []
    DSs = []
    xs = []
    ys = []

    for x, y in zip(cntrs['x'], cntrs['y']):
        conf_sing = SingleConfig(30, rf1, rf2, x, x, y, y)
        FY = RunSingleTactic().run(conf_sing).failure_year
        FYs.append(FY)
        DSs.append(x+y)
        xs.append(x)
        ys.append(y)
    
    out = dict(FY=FYs, DS=DSs, x=xs, y=ys)
    return out




if __name__=="__main__":

    LOAD = True

    rf1, rf2 = 1e-7, 1e-3

    n_doses = 51
    # n_doses = 6

    primary_inoc = dict(RR=rf1*rf2, RS=rf1, SR=rf2, SS=1-rf1-rf2-rf1*rf2)

    
    conf_grid = GridConfig(30, rf1, rf2, n_doses)
    conf_grid.load_saved = LOAD
    conf_grid.primary_inoculum = primary_inoc
    conf_grid.add_string()


    output = RunGrid().run(conf_grid)

    pars = PretendPars(primary_inoc)

    n_conts = 100
    # n_conts = 10

    contour_data = get_contour_data(n_conts, rf1, rf2, load_saved=LOAD)

    DoseSpaceOverview(output, contour_data, conf_grid.config_string_img)
