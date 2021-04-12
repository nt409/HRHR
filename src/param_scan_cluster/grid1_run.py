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

from .functions import ParamScanGrid
from .config import config


def main(config, index):
    ParamScanGrid(config).run(index)


if __name__=="__main__":
    
    if len(sys.argv)==2:
        index = int(sys.argv[1])
        main(config, index)