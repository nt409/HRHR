# proper/fast/prev_run
run_type = "proper"

baseline_dec = 1.11*10**(-2)

baseline = {
    "RF_single": [-10, -3],
    "RF_double": [-15, -4],

    "omega": [0, 1],

    "theta": [0, 12],

    "decay_rate": [(1/3)*baseline_dec, 3*baseline_dec],

    "SR": [0, 1],

    "load_saved": False,
    "save": False,
    "folder_save": "./param_scan/outputs",
}


if run_type == "proper":
    run_pars = {
        "grid_number": 51,
        "n_cont_points": 100,
        "n_years": 45,
        "n_iterations": 5,
    }

elif run_type == "fast":
    run_pars = {
        "grid_number": 5,
        "n_cont_points": 3,
        "n_years": 45,
        "n_iterations": 2,
    }

elif run_type == "prev_run":
    run_pars = {
        "grid_number": 51,
        "n_cont_points": 5,
        "n_years": 35,
        "n_iterations": 5,
    }


def get_par_str(config):
    string = ""

    for key in config.keys():
        attribute = config[key]

        if key in ["load_saved", "save", "folder_save"]:
            continue

        if type(attribute) == float or type(attribute) == int:
            if "n_" in str(key)[:3]:
                key_str = str(key)[2:6]
            else:
                key_str = str(key)[:3]

            string += f"{key_str}={attribute}_"
        else:
            string += f"{str(key)[:3]}_bd={attribute[0]},{attribute[1]}_"

    string = string.replace("0.", "")
    string = string.replace(".", ",")
    string = string.replace("__", "_")

    if string[-1] == "_":
        string = string[:-1]

    return string


config_rand = {**baseline, **run_pars}
config_rand['par_str'] = get_par_str(config_rand)
