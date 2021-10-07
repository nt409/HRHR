import pandas as pd

from .configs import config_app


def combine_haploid_df_app(n_its, n_variants, n_sex_props, n_doses, double_freq_factor_lowest):

    df = pd.DataFrame()

    for index in range(n_its):
        filename = f"./sr_hap/outputs/single/df_app_{n_variants}_{n_sex_props}_{n_doses}_{double_freq_factor_lowest}_{index}.csv"
        tmp = pd.read_csv(filename)
        df = pd.concat([df,tmp])

    filename = f"./sr_hap/outputs/combined/df_app_{n_variants}_{n_sex_props}_{n_doses}_{double_freq_factor_lowest}.csv"
    print(f"Saving df to: {filename}")
    df.to_csv(filename, index=False)

    return df



if __name__=="__main__":

    df = combine_haploid_df_app(**config_app)



