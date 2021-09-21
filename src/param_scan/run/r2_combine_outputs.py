"""Combine the dataframes."""

import pandas as pd

from param_scan.fns.config import config_rand



def main(config, seeds, output_type):

    df = pd.DataFrame()
    
    par_str = config['par_str']

    folder = config['folder_save']

    for seed in seeds:

        temporary = pd.read_csv(f"{folder}/par_scan/{output_type}_seed={seed}_{par_str}.csv")

        df = df.append(temporary, ignore_index=True)

    
    combined_filename = f"{folder}/combined/output_{output_type}_{par_str}.csv"
    print(f"saving combined output to {combined_filename}")
    df.to_csv(combined_filename)





if __name__=="__main__":
    # N_ITS = config_rand['n_iterations']
    N_ITS = 100

    seeds = list(range(N_ITS))
    
    # output_type = "par_df"
    # output_type = "summary_df"
    # main(config_rand, seeds, output_type)

    main(config_rand, seeds, "summary_df")
    main(config_rand, seeds, "par_df")