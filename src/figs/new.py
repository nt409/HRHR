from model.simulator import RunGrid, RunSingleTactic
from model.config_classes import GridConfig, SingleConfig

n_years = 30
n_doses = 11

config = GridConfig(
    n_years,
    None,
    None,
    n_doses,
    primary_inoculum=dict(
        RR=1e-5, RS=1e-3,
        SR=1e-5, SS=1-1e-5-1e-3-1e-5)
)

fungicide_params = dict(
    theta_1=9,
    theta_2=9,
    delta_1=1.11e-2,
    delta_2=1.11e-2,
    omega_1=0.48,
    omega_2=1,
)

model_output = RunGrid(fungicide_params).run(config)

print(model_output)
