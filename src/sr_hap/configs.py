# n_sex_props = 11
# n_doses = 21

n_sex_props = 2
n_doses = 3



config_app = dict(
    n_variants = 3,
    n_its = 27,
    double_freq_factor_lowest = 1e-4,
    
    n_sex_props = n_sex_props,
    n_doses = n_doses,
    )


config_res = dict(
    n_its = 1,
    double_freq_factors = [1e-5, 1, 1e5],

    n_sex_props = n_sex_props,
    n_doses = n_doses,
    )