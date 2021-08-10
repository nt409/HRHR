"""
Combine the dataframes.
"""
from ..fns.config import config_rand
from ..fns.post_process import combine_PS_rand_outputs



def main(config, seeds):
    combine_PS_rand_outputs(config, seeds)


if __name__=="__main__":
    seeds = list(range(32))
    # seeds = [0]
    main(config_rand, seeds)