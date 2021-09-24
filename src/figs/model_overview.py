from model.simulator import RunSingleTactic
from model.config_classes import SingleConfig
from plotting.paper_figs import CombinedModelPlot


# config_sing = SingleConfig(10, 10**(-3), 10**(-6), 1, 1, 0.5, 0.5)

config_sing = SingleConfig(10, None, None, 1, 1, 0.5, 0.5)
config_sing.load_saved = False
RR, RS, SR = (10**(-8), 10**(-3), 10**(-6))

config_sing.primary_inoculum = dict(RR = RR,
    RS = RS,
    SR = SR,
    SS = 1 - RR - RS - SR)
config_sing.add_string()

run_s = RunSingleTactic()
run_s.yield_stopper = 0
output = run_s.run(config_sing)
CombinedModelPlot(output, config_sing.config_string_img)

