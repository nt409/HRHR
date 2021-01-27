from .config_classes import SingleConfig, GridConfig

n_years = 15

res_f1 = 10**(-1)
res_f2 = 10**(-7)

d11 = 1
d12 = 1
d21 = 1
d22 = 1

n_grid = 3


ConfigSingleRun = SingleConfig(n_years, res_f1, res_f2, d11, d12, d21, d22)
ConfigGridRun    =  GridConfig(n_years, res_f1, res_f2, n_grid)

print(ConfigSingleRun.config_string)