"""
Scan over: 
- IRFs
- fung asymptote
- fung decay rate
- partial res
- levels of SR

Process
- run on cluster,
- one seed for each chunk, 
- compile into one df using combine_outputs.py

"""

import sys

from ..fns.calc_full_scan import ParameterScan
from ..fns.config import config_rand


def main(config, seed):
    ParameterScan(config, seed).run()


if __name__=="__main__":
    
    if len(sys.argv)==2:
        seed = int(sys.argv[1])
        main(config_rand, seed)
    else:
        raise Exception("Supply one argument: a random seed")