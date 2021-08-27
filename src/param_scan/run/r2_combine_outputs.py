"""Combine the dataframes."""

from param_scan.fns.config import config_rand
from param_scan.fns.post_process import combine_PS_rand_outputs



def main(config, seeds):
    combine_PS_rand_outputs(config, seeds)


if __name__=="__main__":
    seeds = list(range(32))
    main(config_rand, seeds)