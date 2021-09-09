from model.params import PARAMS
from model.simulator import RunGrid
from model.config_classes import GridConfig
from plotting.paper_figs import DoseSpaceScenariosPlot




def get_data(load_saved, n_doses):

    rf1s, rf2s = 1e-7, 1e-7
    primary_inoc_same = dict(RR=rf1s*rf2s, RS=rf1s, SR=rf2s, SS=1-rf1s-rf2s-rf1s*rf2s)

    rf1d, rf2d = 1e-7, 1e-3
    primary_inoc_diff = dict(RR=rf1d*rf2d, RS=rf1d, SR=rf2d, SS=1-rf1d-rf2d-rf1d*rf2d)



    # same RFs, same fung
    conf_grid = GridConfig(30, None, None, n_doses)
    conf_grid.load_saved = load_saved
    conf_grid.primary_inoculum = primary_inoc_same
    conf_grid.add_string()

    output_SS = RunGrid().run(conf_grid)



    # diff RFs, same fung
    conf_grid = GridConfig(30, None, None, n_doses)
    conf_grid.load_saved = load_saved
    conf_grid.primary_inoculum = primary_inoc_diff
    conf_grid.add_string()

    output_DS = RunGrid().run(conf_grid)






    # same RFs, diff fung
    conf_grid = GridConfig(30, None, None, n_doses)
    conf_grid.load_saved = load_saved
    conf_grid.primary_inoculum = primary_inoc_same
    conf_grid.bs_sex_prop = 1
    conf_grid.add_string()

    fcide_pars = dict(omega_1=1,
                    omega_2=0.95,
                    theta_1=PARAMS.theta_1,
                    theta_2=7.5,
                    delta_1=PARAMS.delta_1,
                    delta_2=PARAMS.delta_2)

    fung_par_str = (f"_fung_pars="
                f"{fcide_pars['omega_1']},{fcide_pars['omega_2']}," 
                f"{fcide_pars['theta_1']},{fcide_pars['theta_2']}")

    conf_grid.config_string = (conf_grid.config_string.split(".pickle")[0] +
                                f"{fung_par_str.replace('.', ',')}.pickle")
    
    output_SD = RunGrid(fcide_pars).run(conf_grid)




    # diff RFs, diff fung
    conf_grid = GridConfig(30, None, None, n_doses)
    conf_grid.load_saved = load_saved
    conf_grid.primary_inoculum = primary_inoc_diff
    conf_grid.add_string()

    fcide_pars = dict(omega_1=1,
                    omega_2=0.8,
                    theta_1=PARAMS.theta_2,
                    theta_2=7,
                    delta_1=PARAMS.delta_1,
                    delta_2=PARAMS.delta_2)

    fung_par_str = (f"_fung_pars="
                f"{fcide_pars['omega_1']},{fcide_pars['omega_2']}," 
                f"{fcide_pars['theta_1']},{fcide_pars['theta_2']}")
    
    conf_grid.config_string = (conf_grid.config_string.split(".pickle")[0] +
                                f"{fung_par_str.replace('.', ',')}.pickle")

    output_DD = RunGrid(fcide_pars).run(conf_grid)

    return output_SS, output_SD, output_DS, output_DD, conf_grid.config_string_img





    
if __name__=="__main__":

    load_saved = True
    n_doses = 51

    data_SS, data_SD, data_DS, data_DD, conf_str_img = get_data(load_saved=load_saved, n_doses=n_doses)

    DoseSpaceScenariosPlot(data_SS, data_SD, data_DS, data_DD, conf_str_img)

