# proper/fast/prev_run
run_type = "fast"

baseline_dec = 1.11*10**(-2)

baseline = {
        "RFS1": [-8,-2],
        "RFS2": [-8,-2],
        "RFD": [-15,-3],

        "asym1": [0.4, 1],
        "asym2": [0.4, 1],
        
        "theta": 12,

        "dec_rate1": [(1/3)*baseline_dec, 3*baseline_dec],
        "dec_rate2": [(1/3)*baseline_dec, 3*baseline_dec],

        "SR": [0,1],
        
        "load_saved": False,
        "save": False,
        "folder_save": "./param_scan/outputs",
        }


if run_type == "proper":
    run_pars = {
        "grid_number": 51,
        "n_cont_points": 71,
        "n_years": 40,
        "NIts": 15,
        }

elif run_type == "fast":
    run_pars = {
        "grid_number": 5,
        "n_cont_points": 3,
        "n_years": 35,
        "NIts": 2,
        }

elif run_type == "prev_run":
    run_pars = {
        "grid_number": 51,
        "n_cont_points": 5,
        "n_years": 35,
        "NIts": 5,   
        }




def get_par_str(config):
    string = ""
    
    for key in config.keys():
        attribute = config[key]

        if key in ["load_saved", "save", "folder_save"]:
            continue

        if type(attribute)==float or type(attribute)==int:
            string += f"{str(key)[:2]}={attribute}"
        else:
            string += f"{str(key)[:2]}=L{attribute[0]},U{attribute[1]}"

    string = string.replace(".", ",")

    return string




config_rand = {**baseline, **run_pars}
config_rand['par_str'] = get_par_str(config_rand)