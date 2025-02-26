"""
Get table with high vs low output
"""

from param_scan.fns.config import config_rand
from param_scan.fns.post_process import PostProcess


def main(config):

    PP = PostProcess(config['par_str'])
    PP.check_high_or_low_dose()
    # do something


if __name__ == "__main__":
    main(config_rand)
