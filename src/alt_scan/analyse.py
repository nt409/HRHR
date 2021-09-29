import itertools
import pandas as pd

from alt_scan.utils import get_alt_scan_params




def main(n_doses, n_its):
    df = pd.DataFrame()

    for alt_strat, index in itertools.product(["alt_12", "alt_21"], list(range(n_its**3))):

        _, _, ii, jj, kk = get_alt_scan_params(n_its, index)
        print(alt_strat, ii, jj, kk)

        filename = f"./alt_scan/outputs/single/out_{n_its}_{n_doses}_{False}_{alt_strat}_{ii}_{jj}_{kk}.csv"
        
        new_df = pd.read_csv(filename)
        new_df["alt_strat"] = alt_strat

        df = pd.concat([df, new_df])



    # print(df.head(-10))
    # print(df.describe())

    print(df[df["asex_alt"]>df["asex_mix"]])
    print(df[df["sex_alt"]>df["sex_mix"]])


    filename = f"./alt_scan/outputs/combined/out_{n_its}_{n_doses}_{alt_strat}.csv"
    print(f"saving df to: {filename}")
    df.to_csv(filename, index=False)





if __name__=="__main__":
    
    n_doses = 21
    n_its = 5

    main(n_doses, n_its)


