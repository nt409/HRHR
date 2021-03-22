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

import copy
import sys

from .functions import ParamScan, get_par_str
from .config import config


def main(config, index):

    config_use = copy.copy(config)

    config_use["SR"] = [config["SR"][index]]

    df = ParamScan(config_use).run_param_scan()
    
    par_str = get_par_str(config_use)

    df.to_csv(f"../outputs/csvs/param_scan_{par_str}.csv", index=False)


if __name__=="__main__":
    
    if len(sys.argv)==2:
        index = int(sys.argv[1])
        main(config, index)