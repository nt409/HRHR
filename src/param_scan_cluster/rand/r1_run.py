"""
Scan over: 
- IRFs
- fung asymptote
- fung decay rate
- partial res
- levels of SR

We use one SR for each chunk, run on cluster,
then compile into one df using post_process.py

"""

import sys

from ..functions import ParamScanRandRFB, ParamScanRandFYY
from ..config import config_rand


def main(config, seed):
    if config['contour_type']=="RFB":
        ParamScanRandRFB(config).run(seed)
    else:
        ParamScanRandFYY(config).run(seed)


if __name__=="__main__":
    
    if len(sys.argv)==2:
        seed = int(sys.argv[1])
        main(config_rand, seed)
    else:
        raise Exception("Supply one argument: a random seed")