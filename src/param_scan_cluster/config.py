import numpy as np

config_grid = {
"RFS1": [10**(k) for k in [-8,-5]],
"RFS2": [10**(k) for k in np.linspace(-8,-1, 5)],
"RFD": [10**(k) for k in np.linspace(-15,-1, 4)],

"asym1": np.linspace(0.5, 1, 3),
"asym2": np.linspace(0.5, 1, 6),

"dec_rate1": [(1.11*10**(-2))*(2**j) for j in [-2,0,2]],
"dec_rate2": [(1.11*10**(-2))*(2**j) for j in [-2,-1,0,1,2]],

"SR": list(np.linspace(0,1,5)) + [0.9, 0.95, 0.99],
"NDoses": 11,
"load_saved": True
}


config_rand = {
    "RFS1": [-8,-2],
    "RFS2": [-8,-2],
    "RFD": [-15,-3],

    "asym1": [0.4, 1],
    "asym2": [0.4, 1],

    "dec_rate1": [0.5*1.11*10**(-2), 2*1.11*10**(-2)],
    "dec_rate2": [0.5*1.11*10**(-2), 2*1.11*10**(-2)],

    "SR": [0,1],
    
    "grid_number": 26,

    "NIts": 75,
    "load_saved": False,
    
    "contour_type": "RFB",
}

# config_rand = {
#     "RFS1": [-8,-2],
#     "RFS2": [-8,-2],
#     "RFD": [-15,-3],

#     "asym1": [0.4, 1],
#     "asym2": [0.4, 1],

#     "dec_rate1": [0.5*1.11*10**(-2), 2*1.11*10**(-2)],
#     "dec_rate2": [0.5*1.11*10**(-2), 2*1.11*10**(-2)],

#     "SR": [0,1],
    
#     "grid_number": 21,

#     "NIts": 1500,
#     "load_saved": True,
    
#     "contour_type": "RFB",
# }