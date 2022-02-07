"""
Scan over: 
- IRFs
- fung asymptote
- fung curvature
- fung decay rate
- partial res
- levels of between-season SR

Process
- run on cluster
- one seed for each chunk
- compile into one df using r2_combine_outputs.py
"""

import sys

from param_scan.fns.calc_full_scan import ParameterScan
from param_scan.fns.config import config_rand


def main(config, seed):
    ParameterScan(config, seed).run()


if __name__ == "__main__":

    if len(sys.argv) != 2:
        raise Exception("Supply one argument: a random seed")

    seed = int(sys.argv[1])
    main(config_rand, seed)
