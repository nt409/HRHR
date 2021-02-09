from .config_classes import SingleConfig, GridConfig

n_years = 20

res_f1 = 10**(-3)
res_f2 = 10**(-3)

d11 = 0.9
d12 = 0.9
d21 = 0.9
d22 = 0.9

n_grid = 7


ConfigSingleRun = SingleConfig(n_years, res_f1, res_f2, d11, d12, d21, d22)
ConfigGridRun    =  GridConfig(n_years, res_f1, res_f2, n_grid)