import numpy as np
from math import floor
##
import sys
# sys.path.insert(0, '/utils/')
##
from utils.functions_HRHR import master_loop_grid_of_tactics, object_dump, object_open, cluster_chunk
from utils.parameters_HRHR import params
#----------------------------------------------------------------------------------------------
asex_dict = object_open(params.JSON_path+'global.json','json')
##
asex_dict['phi_vec']   = np.linspace(asex_dict['log_phi_min'],0,asex_dict['n_phi'])
asex_dict['phi_vec_rr']   = np.linspace(2*asex_dict['log_phi_min'],0,2*asex_dict['n_phi'])
param_string_rec  = ',n_phi=' + str(asex_dict['n_phi']) + ',n_d=' + str(asex_dict['n_d']) + ',n_rec=' + str(asex_dict['n_recurs'])
#----------------------------------------------------------------------------------------------
if len(sys.argv)==3: # to run write HR_HR/Scripts/run2_optimal_asexual_cluster.py phi_con_0 2 ##or phi_con_x for job number x# index is 0,1,...,floor(N/16); so that N=47 gives 0,1,2, N=48 gives 0,1,2,3
    ##
    x = sys.argv[1]
    y = int(x[8:])
    index = int(sys.argv[2])
    if y<32: # 0 to 31
        for i in [index]:#range(1+floor(2*asex_dict['n_phi']/32)):
            z = y + 32*i # 32*index + [y; number in 0,1,...,31]
            z = int(z)
            if z<2*asex_dict['n_phi']: # z in 0,1,...,2*N_phi - 1
                config_phi_rr = params.JSON_path + 'phi_con/phi_con_' + str(z) + '.json' # so need 2*n_phi to be a multiple of 32
                dict_phi_rr = object_open(config_phi_rr,'json')
                phi_index = [dict_phi_rr['phi_rr']]
                if phi_index[0] <= (2*asex_dict['n_phi'] - 1):
                    cluster_chunk(phi_index[0],asex_dict,param_string_rec)
    