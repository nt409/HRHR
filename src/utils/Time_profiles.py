import timeit
# local setup
import numpy as np
from math import ceil, log, floor, log10, exp
import pickle
import json
import sys
from Functions_and_plotting.functions_HRHR import master_loop_grid_of_tactics
from Functions_and_plotting.parameters_HRHR import params, params_dict
# from Optimal_HRHR_asexual_cluster import cluster_chunk
#----------------------------------------------------------------------------------------------
# with open(params.JSON_path+'global.json') as config_file:
#     asex_dict = json.load(config_file)
# ##
# asex_dict['phi_vec']   = np.linspace(asex_dict['log_phi_min'],0,asex_dict['n_phi'])
# #----------------------------------------------------------------------------------------------
# param_string_rec  = ',n_phi=' + str(asex_dict['n_phi']) + ',n_d=' + str(asex_dict['n_d']) + ',n_rec=' + str(asex_dict['n_recurs'])


h2 = {'hi':5}
h3 = {**h2, **params_dict}
#----------------------------------------------------------------------------------------------
# t = timeit.timeit("params.Fung1_ind']",'from parameters_HRHR import params',number=1000)
# t2 = timeit.timeit("params2.Fung1_ind",'from parameters_HRHR import params2',number=1000)
# t3 = timeit.timeit("params3.Fung1_ind",'from parameters_HRHR import params3',number=1000)
# print(t)
# print(t2)
# print(t3)

prr = 0.1
prs = 0.1
psr = 0.1
pss = 1-prr-prs-psr
n_d=4
t4 = timeit.timeit("master_loop_grid_of_tactics(4,1,p_rr=0.1,p_rs=0.1,p_sr=0.1,p_ss=0.7)",'from calculator_HRHR import master_loop_grid_of_tactics',number=100)
print(t4)