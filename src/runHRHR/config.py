from .config_classes import SingleConfig, GridConfig

n_years = 20

res_f1 = 10**(-5)
res_f2 = 10**(-5)

d11 = 1
d12 = 1
d21 = 0.5
d22 = 0.5


n_grid = 4

ConfigSingleRun = SingleConfig(n_years, res_f1, res_f2, d11, d12, d21, d22)
ConfigGridRun    =  GridConfig(n_years, res_f1, res_f2, n_grid)