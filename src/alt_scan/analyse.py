import pandas as pd

def main(n_doses, n_its):
    filename = f"./alt_scan/outputs/combined/out_{n_its}_{n_doses}.csv"
    df = pd.read_csv(filename)

    print(df.head())
    print(df.describe())

    print(df.loc[df["asex_mix"]<df["sex_mix"], ["sex_mix"]].count())
    print(df.loc[df["asex_mix"]==df["sex_mix"]])





if __name__=="__main__":
    
    n_doses = 21
    n_its = 5

    main(n_doses, n_its)