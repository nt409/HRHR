import pandas as pd


def combine_haploid_df(n_its, n_sex_props, n_doses, double_freq_factors):
    df = pd.DataFrame()

    dff_str = ",".join([str(ee) for ee in double_freq_factors])
    
    for ii in range(3*n_its):
        filename = f"./sr_hap/outputs/single/df_{n_its}_{n_sex_props}_{n_doses}_{dff_str}_{ii}.csv"
        tmp = pd.read_csv(filename)
        df = pd.concat([df,tmp])

    filename = f"./sr_hap/outputs/combined/df_{n_its}_{n_sex_props}_{n_doses}_{dff_str}.csv"
    print(f"Saving df to: {filename}")
    df.to_csv(filename, index=False)

    return df



if __name__=="__main__":
    n_its = 10
    n_sex_ps = 11
    n_doses = 21

    double_freq_factors = [1e-5, 1, 1e5]
    df = combine_haploid_df(n_its, n_sex_ps, n_doses, double_freq_factors)

    print(df)



