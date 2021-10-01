import pandas as pd


def combine_haploid_df(n_variants, n_sex_props, n_doses, double_freq_factor_lowest):
    df = pd.DataFrame()

    for index in range(n_its):
        filename = f"./sr_hap/outputs/single/df_{n_variants}_{n_sex_props}_{n_doses}_{double_freq_factor_lowest}_{index}.csv"
        tmp = pd.read_csv(filename)
        df = pd.concat([df,tmp])

    filename = f"./sr_hap/outputs/combined/df_{n_variants}_{n_sex_props}_{n_doses}_{double_freq_factor_lowest}.csv"
    print(f"Saving df to: {filename}")
    df.to_csv(filename, index=False)

    return df



if __name__=="__main__":
    n_variants = 3
    n_sex_ps = 11
    n_doses = 21
    n_its = 27
    double_freq_factor_lowest = 1e-4

    df = combine_haploid_df(n_its, n_variants, n_sex_ps, n_doses, double_freq_factor_lowest)

    print(df)



