import numpy as np

# config = {
# "RFS1": [10**(k) for k in [-5]],
# "RFS2": [10**(k) for k in [-5]],
# "RFD": [10**(k) for k in [-6, -5]],

# "asym1": [0.8],
# "asym2": [0.8],

# "dec_rate1": [(1.11*10**(-2))*(2**j) for j in [0]],
# "dec_rate2": [(1.11*10**(-2))*(2**j) for j in [0]],

# "SR": [0,0.5,1],
# "NDoses": 3
# }

config = {
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
    "RFS1": [-8,-1],
    "RFS2": [-8,-1],
    "RFD": [-15,-3],

    "asym1": [0.4, 1],
    "asym2": [0.4, 1],

    "dec_rate1": [0.5*1.11*10**(-2), 2*1.11*10**(-2)],
    "dec_rate2": [0.5*1.11*10**(-2), 2*1.11*10**(-2)],

    "SR": [0,1],
    "NDoses": 5,
    "NIts": 100,
    "load_saved": True
}