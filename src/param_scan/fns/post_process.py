from param_scan.fns.calc_method_contour import RunAlongContourDFs
import pandas as pd
import copy
import numpy as np


from model.simulator import RunGrid
from param_scan.fns.pars import RandomPars
from plotting.figures import DoseSpaceScenarioSingle

# TOC
# combine_PS_rand_outputs
# PostProcess
# ProcessedDF



def combine_PS_rand_outputs(config, seeds, output_type):

    df = pd.DataFrame()
    
    par_str = config['par_str']

    folder = config['folder_save']

    for seed in seeds:

        temporary = pd.read_csv(f"{folder}/par_scan/{output_type}_seed={seed}_{par_str}.csv")

        df = df.append(temporary, ignore_index=True)

    
    combined_filename = f"{folder}/combined/output_{output_type}_{par_str}.csv"
    print(f"saving combined output to {combined_filename}")
    df.to_csv(combined_filename)





class PostProcess:

    def __init__(self, par_str):
        self.folder = "./param_scan/outputs"

        df_in = pd.read_csv(f"{self.folder}/combined/output_summary_df_{par_str}.csv")

        self.df = df_in.drop(["Unnamed: 0"], axis=1)
        
        self.par_df = pd.read_csv(f"{self.folder}/combined/output_par_df_{par_str}.csv")

        self.processed_df = ProcessedDF(self.folder, self.par_df, self.df).df

        self.par_str = par_str
        






    


    def analyse_processed_df(self):

        df = copy.copy(self.processed_df)

        out = pd.DataFrame()

        for strat in ["c_R_maxCont%", "c_E_maxCont%", "c_E_lowDoseMax%"]:

            row = self.get_outcome_row_this_strat(df, strat)
            out = out.append(row, ignore_index=True)
        
        return out




    def get_outcome_row_this_strat(self, df, strategy):

        mean = df[strategy].mean()

        conditional_mean_less_than = df[df[strategy]<100][strategy].mean()
        
        conditional_mean_leq = df[df[strategy]<=100][strategy].mean()

        geq = df[df[strategy]>=100].shape[0]
        greater = df[df[strategy]>100].shape[0]
        lesser = df[df[strategy]<100].shape[0]
        equal = df[df[strategy]==100].shape[0]

        total = df.shape[0]

        sum_total = sum([greater, lesser, equal])

        out = dict(
                    geq=geq, 
                    greater=greater,
                    lesser=lesser,
                    equal=equal,
                    total=total,
                    sum_total=sum_total,
                    work_pc=round(100*geq/sum_total,1),
                    mean=round(mean,1),
                    conditional_mean_less_than=round(conditional_mean_less_than,1),
                    conditional_mean_leq=round(conditional_mean_leq,1),
                    strategy=strategy
                    )
        
        return out






    def get_failed_pars(self):
        df = self.processed_df

        failed = df[df['c_R_maxCont%']<100]

        runs_that_failed = failed["run"].unique()
        failed_pars = self.get_params_for_specific_runs(runs_that_failed)

        return failed_pars

    





    def get_params_for_specific_runs(self, which_runs):

        out = pd.DataFrame()
        for rr in which_runs:
            this_run = self.processed_df[self.processed_df["run"]==rr].iloc[0,:]
            out = out.append(this_run, ignore_index=True)
        
        return out








        




    def check_high_or_low_dose(self):

        my_df = copy.copy(self.df)

        my_df['high_better_than_low'] = my_df['c_R_highDoseMaxEL'] >= my_df['c_R_lowDoseMaxEL']

        # strats = ["minEqDose", "fullDose"]
        
        # for string in strats:
        #     my_df[string + "%"] = 100*my_df[string + "EL"]/my_df["max_grid_EL"]
        
        grouped = my_df.groupby(["run"]).first()

        df_out = pd.DataFrame(grouped)

        df_out = df_out.reset_index()
        
        df_out = df_out.sort_values(['high_better_than_low', 'sr_prop'])

        
        sex_and_high_eff = df_out[((df_out['sr_prop']>0.9)
                            & (df_out['omega_1']>0.9) 
                            & (df_out['omega_2']>0.9) 
                            # & (df_out['delta_1']>0.01)
                            # & (df_out['delta_2']>0.01)  
                            # & (df_out['RS']<0.00001)
                            # & (df_out['SR']<0.00001)
                            )]
        
        # sex_and_high_eff = df_out[(df_out['sr_prop']>0.6)]
        
        print("\n")
        print("worked/total:",sex_and_high_eff['high_better_than_low'].sum(), sex_and_high_eff.shape[0])
        print("\n")
        print(sex_and_high_eff[['high_better_than_low', 'sr_prop',
                'omega_1', 'omega_2', 'delta_1', 'delta_2',
                'RR'
                ]].sort_values(['high_better_than_low', 'sr_prop']))


        # print(df_out.loc[~df_out['high_better_than_low']]['sr_prop'].mean())
                
        filename = f"{self.folder}/par_scan/high_or_low_dose_{len(df_out)}.csv"

        print(f"Saving high or low dose csv to: \n{filename}")
        
        df_out.to_csv(filename)




  
    
    
    def re_run_grid(self, NDoses, run_indices):

        df_test = self.get_params_for_specific_runs(run_indices)

        for ii in range(df_test.shape[0]):

            pars = df_test.iloc[int(ii),:]
           
            print("\nRe-running run:", df_test.iloc[int(ii),:].run, "\n")

            rp = self._get_RPs(pars, NDoses)
            
            grid_output = RunGrid(rp.fung_parms).run(rp.grid_conf)

            conf_str = rp.grid_conf.config_string_img

            FY = grid_output.FY
            opt_region = FY == np.amax(FY)
                
            n_opt_doses = opt_region.sum()

            print(f"Number of optimal dose combos: {n_opt_doses}")

            # plot output
            DoseSpaceScenarioSingle(grid_output, conf_str)





    def re_run_cont(self, NDoses, N_cont_doses, run_indices):

        df_test = self.get_params_for_specific_runs(run_indices)

        out = pd.DataFrame()

        for ii in range(df_test.shape[0]):

            pars = df_test.iloc[int(ii),:]

            this_run_ind = int(df_test.iloc[int(ii),:].run)

            print(f"\nRe-running run: {this_run_ind} \n")

            rp = self._get_RPs(pars, NDoses)

            grid_output = RunGrid(rp.fung_parms).run(rp.grid_conf)
           
            RFB_dfs = RunAlongContourDFs(rp, grid_output, N_cont_doses, "RFB")


            data = dict(run=this_run_ind,
                    c_R_maxContEL = RFB_dfs.summary.c_R_maxContEL,
                    max_grid_EL = np.amax(grid_output.FY)
                    )

            out = out.append(data, ignore_index=True)
            
            # plot output
            # conf_str = rp.grid_conf.config_string_img
            # DoseSpaceScenarioSingle(grid_output, conf_str)


        for col in ["c_R_maxContEL", "max_grid_EL"]:
            out[col] = out[col].astype("int")

        out["worked"] = out["c_R_maxContEL"] >= out["max_grid_EL"]


        ind_str = ",".join([str(rr) for rr in run_indices])
        filename = f"{self.folder}/par_scan/re_run_{NDoses}_{N_cont_doses}_{ind_str}.csv"
        print(f"Saving re-run to: {filename}")
        print(out)
        out.to_csv(filename, index=False)





    @staticmethod
    def _get_RPs(pars, NDoses):

        config = {'load_saved': True, 'save': True, 'n_years': 35}

        RP = RandomPars(config, None)

        RP.get_inoc_dict(pars["RS"], pars["SR"], pars["RR"])
        
        RP.fung_parms = RP.get_fung_parms_dict(pars["omega_1"], pars["omega_2"], 
                                pars["delta_1"], pars["delta_2"], pars["theta_1"])
        
        RP.sr_prop = pars["sr_prop"]

        RP.path_and_fung_pars = (pars["RS"], pars["SR"], pars["RR"], pars["omega_1"],
                    pars["omega_2"], pars["delta_1"], pars["delta_2"])

        RP.get_grid_conf(NDoses)

        return RP














class ProcessedDF:
    """
    Performs some basic operations on the data and returns a dataframe with
    parameters and interesting output info only.
    """

    def __init__(self, folder, par_df, df_input):
        self.folder = folder
        
        df_inter = self._process_df(par_df, df_input)

        self.df = self._check_if_never_failed(df_inter)
        
        # self._save_df()




    def _process_df(self, par_df, data_in):

        par_df.drop(["Unnamed: 0"], axis=1, inplace=True)

        good_cols = [
                    "run",
                    "c_R_maxContEL",
                    "c_E_maxContEL",
                    "c_E_lowDoseMaxEL",
                    "c_R_lowDoseMaxEL",
                    "c_R_medDoseMaxEL",
                    "c_R_highDoseMaxEL",
                    "max_grid_EL",
                    "O_corner_01",
                    "O_corner_10",
                    "O_fullDoseEL",
                    "I_R_best_value",
                    "I_E_best_value",
                    "c_R_min_opt_dist_from_contour",
                    "c_E_min_opt_dist_from_contour",
                    ]

        data_use = data_in.loc[:, good_cols]


        left = par_df.set_index(["run"], drop=False)
        right = data_use.set_index(["run"])

        data = left.join(right)
        
        
        data['maxAlongContour'] = data['c_R_maxContEL'] >= data['max_grid_EL']

        strats = ["c_R_maxCont", "c_E_maxCont", "c_E_lowDoseMax"]
        
        for string in strats:
            data.loc[:, string + "%"] = 100*data[string + "EL"]/data["max_grid_EL"]
            
        data['min_corner'] = data[["O_corner_01", "O_corner_10"]].min(axis=1)

        data['ERFB_diff_from_opt'] = data.apply(self.get_diff_from_opt_RFB, axis=1)

        data['ESFY_diff_from_opt'] = data.apply(self.get_diff_from_opt_EqSel, axis=1)

        integer_cols = [
                        'run',
                        'min_corner',
                        'ERFB_diff_from_opt',
                        'ESFY_diff_from_opt',
                        'max_grid_EL',
                        'c_R_maxContEL',
                        'I_R_best_value',
                        'c_E_maxContEL',
                        'I_E_best_value',
                        'c_E_lowDoseMaxEL',
                        ]

        for col in integer_cols:
            data[col] = data[col].astype("int")

        return data


    def get_diff_from_opt_RFB(self, data):
        return data['max_grid_EL'] - max(data['c_R_maxContEL'], data['I_R_best_value'])
    
    def get_diff_from_opt_EqSel(self, data):
        return data['max_grid_EL'] - max(data['c_E_maxContEL'], data['I_E_best_value'])
    



    def _check_if_never_failed(self, df):

        # avoid the "-1" case where has never failed:
        if df[df["O_fullDoseEL"]<=0].shape[0]:
            n_never_fail = df[df["O_fullDoseEL"]<=0].shape[0]
            print(f"{n_never_fail} runs never failed - need longer n_years. Filtering them out for now.")
            df = df[df["O_fullDoseEL"]>0]
        
        return df





    def _save_df(self):

        df = self.df
        
        filename = f"{self.folder}/par_scan/processed_{len(df)}.csv"

        print(f"Saving processed data to: \n{filename}")

        df.to_csv(filename)



