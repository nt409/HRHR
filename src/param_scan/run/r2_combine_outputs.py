"""Combine the dataframes."""

from param_scan.fns.config import config_rand
from param_scan.fns.post_process import combine_PS_rand_outputs



def main(config, seeds, output_type):
    combine_PS_rand_outputs(config, seeds, output_type)


if __name__=="__main__":
    seeds = list(range(32))
    output_type = "par_df"
    # output_type = "summary_df"

    main(config_rand, seeds, output_type)