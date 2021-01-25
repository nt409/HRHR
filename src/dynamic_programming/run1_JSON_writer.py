import json
##
# import sys
# sys.path.insert(0, '/utils/')
##
from utils.parameters_HRHR import params
from utils.functions_HRHR import object_dump, object_open, calculator_d_free
#----------------------------------------------------------------------------------------------
asex_dict = object_open(params.JSON_path+'global.json','json')
for i in range(2*asex_dict['n_phi']):
    config_phi_rr = params.JSON_path + 'phi_con/' + 'phi_con_' + str(i) + '.json'
    config = dict(phi_rr=i)
    object_dump(config_phi_rr, config, 'json')