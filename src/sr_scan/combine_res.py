import pandas as pd

from .configs import config_res

def combine_haploid_df_res(n_its, n_sex_props, n_doses, double_freq_factors):

    df = pd.DataFrame()

    dff_str = ",".join([str(ee) for ee in double_freq_factors])

    for index in range(3*n_its):
        filename = f"./sr_scan/outputs/single/df_res_{n_its}_{n_sex_props}_{n_doses}_{dff_str}_{index}.csv"
        tmp = pd.read_csv(filename)
        df = pd.concat([df,tmp])

    filename = f"./sr_scan/outputs/combined/df_res_{n_its}_{n_sex_props}_{n_doses}_{dff_str}.csv"
    print(f"Saving df to: {filename}")
    df.to_csv(filename, index=False)

    return df



if __name__=="__main__":

    df = combine_haploid_df_res(**config_res)

