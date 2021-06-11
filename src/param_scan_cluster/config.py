

# proper/fast/prev_run
run_type = "proper"


if run_type == "proper":
    config_rand = {
        "RFS1": [-8,-2],
        "RFS2": [-8,-2],
        "RFD": [-15,-3],

        "asym1": [0.4, 1],
        "asym2": [0.4, 1],

        "dec_rate1": [0.5*1.11*10**(-2), 2*1.11*10**(-2)],
        "dec_rate2": [0.5*1.11*10**(-2), 2*1.11*10**(-2)],

        "SR": [0,1],
        
        "grid_number": 50,

        "n_cont_points": 50,

        "n_years": 50,

        "NIts": 25,
        
        "load_saved": True,   
    }


elif run_type == "fast":
    config_rand = {
        "RFS1": [-8,-2],
        "RFS2": [-8,-2],
        "RFD": [-15,-3],

        "asym1": [0.4, 1],
        "asym2": [0.4, 1],

        "dec_rate1": [0.5*1.11*10**(-2), 2*1.11*10**(-2)],
        "dec_rate2": [0.5*1.11*10**(-2), 2*1.11*10**(-2)],

        "SR": [0,1],
        
        "grid_number": 11,

        "n_cont_points": 5,

        "n_years": 35,

        "NIts": 5,
        
        "load_saved": True,
    }


elif run_type == "prev_run":
    config_rand = {
        "RFS1": [-8,-2],
        "RFS2": [-8,-2],
        "RFD": [-15,-3],

        "asym1": [0.4, 1],
        "asym2": [0.4, 1],

        "dec_rate1": [0.5*1.11*10**(-2), 2*1.11*10**(-2)],
        "dec_rate2": [0.5*1.11*10**(-2), 2*1.11*10**(-2)],

        "SR": [0,1],
        
        "grid_number": 51,

        # "n_cont_points": 41,

        # "n_years": 35,

        "NIts": 32,
        "load_saved": True,
        
    }