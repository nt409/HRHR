import itertools
import pandas as pd

# from alt_scan.utils import get_alt_scan_params




def main(n_doses, n_its):
    df = pd.DataFrame()

    for alt_strat, index in itertools.product(["alt_12", "alt_21"], list(range(n_its))):

        filename = f"./alt_scan/outputs/single/out_rand_{n_doses}_{alt_strat}_{index}.csv"
        
        new_df = pd.read_csv(filename)
        new_df["alt_strat"] = alt_strat

        df = pd.concat([df, new_df])



    # print(df.head(-10))
    # print(df.describe())

    print("WORKED RUNS; MIX AT LEAST AS GOOD")
    print(df[df["asex_alt"]<=df["asex_mix"]])
    print(df[df["sex_alt"] <= df["sex_mix"]])

    print("FAILED RUNS; ALT OUTPERFORMS MIX")
    print(df[df["asex_alt"]>df["asex_mix"]])
    print(df[df["sex_alt"]>df["sex_mix"]])
    
    print("NB 100 iterations but over alt_21 and alt_12")


    filename = f"./alt_scan/outputs/combined/out_{n_its}_{n_doses}_{alt_strat}.csv"
    print(f"saving df to: {filename}")
    df.to_csv(filename, index=False)





if __name__=="__main__":
    
    n_doses = 21
    n_its = 100

    main(n_doses, n_its)


