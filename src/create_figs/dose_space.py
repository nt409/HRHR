from model.utils import object_dump
import os
import pickle
from model.strategy_arrays import EqualResFreqBreakdownArray, EqualSelectionArray
from model.simulator import RunGrid, RunSingleTactic
from model.config_classes import GridConfig, SingleConfig
from plotting.paper_figs import DoseSpaceScenariosPlot




class PretendPars:
    fung_parms = None
    
    def get_single_conf(self, dose1, dose2):

        conf = SingleConfig(30, 10**(-7), 10**(-3),
                                dose1, dose1, dose2, dose2,
                                primary_inoculum=None)
        
        config_out = self._process_conf(conf)

        self.sing_conf = config_out

        return config_out
    

    def _process_conf(self, conf):

        conf.bs_sex_prop = 0

        conf.load_saved = False
        
        conf.save = False

        conf.add_string()

        # config_out = self._update_par_scan_conf_str(conf)

        return conf





def get_this_contour_data(z, n_contours, rf1, rf2):
    cntrs = z.find_contours(pars, n_contours)

    FYs = []
    DSs = []

    for x, y in zip(cntrs['x'], cntrs['y']):
        conf_sing = SingleConfig(30, rf1, rf2, x, x, y, y)
        FY = RunSingleTactic().run(conf_sing).failure_year
        FYs.append(FY)
        DSs.append(x+y)
    
    out = dict(FY=FYs, DS=DSs)
    return out



def get_contour_data(n_contours, rf1, rf2):
    
    filename = f'../outputs/saved_runs/contours_{n_contours}_{rf1}_{rf2}.pickle'

    if os.path.isfile(filename):
    # if False:
        with open(filename, 'rb') as f:
            loaded = pickle.load(f)
        return loaded
    
    d1 = get_this_contour_data(EqualResFreqBreakdownArray(output), n_contours, rf1, rf2)
    d2 = get_this_contour_data(EqualSelectionArray(output), n_contours, rf1, rf2)
    run = dict(ERFB=d1, ESFY=d2)
    object_dump(filename, run)

    return run



if __name__=="__main__":

    rf1, rf2 = 10**(-7), 10**(-3)
    
    conf_grid = GridConfig(30, rf1, rf2, 51)
    # conf_grid = GridConfig(30, 10**(-7), 10**(-3), 6)

    conf_grid.load_saved = True
    output = RunGrid().run(conf_grid)

    pars = PretendPars()

    contour_data = get_contour_data(100, rf1, rf2)

    DoseSpaceScenariosPlot(output, contour_data, conf_grid.config_string_img)
