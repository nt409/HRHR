import pandas as pd
import copy
import numpy as np


from utils.functions import RunGrid
from plotting.figures import DoseSpaceScenariosPlot
from .pars import RandomPars

# TOC
# combine_PS_rand_outputs
# PostProcess
# MaxAlongContourDF



def combine_PS_rand_outputs(config, seeds):

    df = pd.DataFrame()
    
    par_str = config['par_str']

    folder = config['folder_save']

    for seed in seeds:

        temporary = pd.read_csv(f"{folder}/par_scan/summary_seed={seed}_{par_str}.csv")

        df = df.append(temporary, ignore_index=True)

    df.to_csv(f"{folder}/combined/output_summary_{par_str}.csv")





class PostProcess:

    def __init__(self, folder, par_str):
        df_in = pd.read_csv(f"{folder}/combined/output_summary_{par_str}.csv")
        self.df = df_in.drop(["Unnamed: 0"], axis=1)

        self.folder = folder
        self.par_str = par_str
        







    def get_maximum_along_contour_df(self):
        self.max_along_contour_df = MaxAlongContourDF(self.folder, self.df).df


    











    def analyse_max_contour_df(self):

        df = copy.copy(self.max_along_contour_df)

        # df = df[df['min_corner']>0]

        eq_sel_df = self._get_non_null_df(df, "EqSel_worked_geq")

        self.filtered_dataframe_outcome(df, "RFB_maxCont%")
        self.filtered_dataframe_outcome(df, "fullDose%")
        self.filtered_dataframe_outcome(df, "minEqDose%")

        self.filtered_dataframe_outcome(eq_sel_df, "EqSel_maxCont%")
        self.filtered_dataframe_outcome(eq_sel_df, "EqSel_lowDoseMax%")

        self.check_IVT_method(df, "RFB")
        self.check_IVT_method(eq_sel_df, "EqSel")





    def _get_non_null_df(self, df_in, col_to_check_if_null):
        invalid_runs = df_in[col_to_check_if_null].isnull()
        return df_in[~invalid_runs]



    def filtered_dataframe_outcome(self, df, strategy):

        mean = df[strategy].mean()

        conditional_mean = df[df[strategy]<100][strategy].mean()

        worked = df[df[strategy]>=100].shape[0]
        greater = df[df[strategy]>100].shape[0]
        lesser = df[df[strategy]<100].shape[0]
        equal = df[df[strategy]==100].shape[0]

        total = df.shape[0]

        sum_total = sum([greater, lesser, equal])

        out = dict(worked=worked, 
                    greater=greater,
                    lesser=lesser,
                    equal=equal,
                    total=total,
                    sum_total=sum_total,
                    work_pc=round(100*worked/sum_total,1),
                    mean=round(mean,1),
                    conditional_mean=round(conditional_mean,1),
                    strategy=strategy
                    )
        
        print(out)

        return out



    def check_IVT_method(self, df, method):

        # if False: # method=="EqSel":
        #     IVT_true = df[f'best_region_{method}']=="True"
        #     IVT_false = df[f'best_region_{method}']=="False"
        #     IVT_NA = df[f'best_region_{method}'].isin(["True", "False"])
        # else:
        
        IVT_true = df[f'best_region_{method}']==True
        IVT_false = df[f'best_region_{method}']==False
        IVT_NA = df[f'best_region_{method}'].isin([True, False])

        opt_true = (df[f'{method}_worked_geq']==True)
        opt_false = (df[f'{method}_worked_geq']==False)
        opt_NA = (df[f'{method}_worked_geq'].isnull())

        out = dict(
            IVT_method_T = df[IVT_true].shape[0],
            IVT_method_F = df[IVT_false].shape[0],
            IVT_method_NA = df[~IVT_NA].shape[0],
            
            opt_method_T = df[opt_true].shape[0],
            opt_method_F = df[opt_false].shape[0],
            opt_method_NA = df[opt_NA].shape[0],

            both_succeeded = df[(opt_true & IVT_true)].shape[0],
            both_failed = df[(~opt_true & ~IVT_true)].shape[0],
            
            either_succeeded = df[(opt_true | IVT_true)].shape[0],
            opt_but_not_IVT = df[(opt_true & ~IVT_true)].shape[0],
            IVT_but_not_opt = df[(~opt_true & IVT_true)].shape[0],
            method = method,
            )

        failed = df[((~opt_true) & (~IVT_true))]
        
        failed_runs = failed[[f'max_grid_EL', f'{method}_maxContEL', f'best_value_{method}', 
                    f'{method}_diff_from_opt', 'run']].sort_values(by=[f'{method}_diff_from_opt', 'run'])

        print("\n")
        print(out)
        print("\n")
        print(f"Testing {method}; these runs failed on both methods:")
        print("\n")        
        
        print(failed_runs)






    def analyse_failed(self):

        df = copy.copy(self.max_along_contour_df)

        fail = self._get_failed_runs(df)

        print("\n")
        print("These runs failed:\n")

        print(fail[['RFB_diff_from_opt',
                    'run',
                    'max_grid_EL',
                    'RFB_maxContEL',
                    'best_value_RFB',
                    ]].to_string())

        
    @staticmethod
    def _get_failed_runs(df):
        return df[df['RFB_maxCont%']<100]

    



    def which_runs_worked_max_cont(self):
        df = copy.copy(self.max_along_contour_df)

        failed = df[df['RFB_maxCont%']<100]

        runs_that_failed = failed["run"].unique()

        failed_pars = self.get_params_for_specific_runs(runs_that_failed)

        n_fail = failed_pars.shape[0]

        failed_pars.to_csv(f"{self.folder}/par_scan/failed_{n_fail}.csv")

    





    def get_params_for_specific_runs(self, which_runs):

        par_df = copy.copy(self.df)

        out = pd.DataFrame()
        for rr in which_runs:
            this_run = par_df[par_df["run"]==rr].iloc[0,:]
            out = out.append(this_run, ignore_index=True)
        
        return out








        




    def check_high_or_low_dose(self):

        my_df = copy.copy(self.df)

        my_df['high_better_than_low'] = my_df['RFB_highDoseMaxEL'] >= my_df['RFB_lowDoseMaxEL']

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




  
    
    
    def re_run(self, NDoses, run_indices):

        df_test = self.get_params_for_specific_runs(run_indices)

        for ii in range(df_test.shape[0]):

            pars = df_test.iloc[int(ii),:]
           
            print("\nRe-running run:", df_test.iloc[int(ii),:].run, "\n")

            grid_config, fung_params = self._get_grid_config_and_fung_pars(pars, NDoses)

            grid_output = RunGrid(fung_params).run(grid_config)

            conf_str = grid_config.config_string_img

            FY = grid_output['FY']
            opt_region = FY == np.amax(FY)
                
            n_opt_doses = opt_region.sum()

            print(f"Number of optimal dose combos: {n_opt_doses}")

            # plot output
            # dose_grid_heatmap(grid_output, grid_config, "FY", conf_str)
            
            # eq_RFB_contours(grid_output, grid_config, title=f"Run={str(pars.run)}")
            DoseSpaceScenariosPlot(grid_output, conf_str)






    @staticmethod
    def _get_grid_config_and_fung_pars(pars, NDoses):

        config = {'load_saved': True, 'n_years': 35}

        RP = RandomPars(config, None)

        RP.get_inoc_dict(pars["RS"], pars["SR"], pars["RR"])
        
        fung_params = RP.get_fung_parms_dict(pars["omega_1"], pars["omega_2"], 
                                pars["delta_1"], pars["delta_2"])
        
        RP.sr_prop = pars["sr_prop"]

        RP.path_and_fung_pars = (pars["RS"], pars["SR"], pars["RR"], pars["omega_1"],
                    pars["omega_2"], pars["delta_1"], pars["delta_2"])

        grid_config = RP.get_grid_conf(NDoses)

        return grid_config, fung_params














class MaxAlongContourDF:

    def __init__(self, folder, df_input):
        self.folder = folder
        self.get_and_save(df_input)
    


    def get_and_save(self, df_input):
        df_inter = self._get_intermediate_df(df_input)

        self.df = self._tidy_df(df_inter)
        
        self._save_df()





    def _get_intermediate_df(self, data):

        data.fillna(0)

        data = pd.DataFrame(data)
        
        data['maxAlongContour'] = data['RFB_maxContEL'] >= data['max_grid_EL']

        strats = ["RFB_maxCont", "EqSel_maxCont", "minEqDose", "fullDose", "EqSel_lowDoseMax"]
        
        for string in strats:
            data.loc[:, string + "%"] = 100*data[string + "EL"]/data["max_grid_EL"]
            
        data['min_corner'] = data[["corner_01", "corner_10"]].min(axis=1)

        data['RFB_diff_from_opt'] = data.apply(self.get_diff_from_opt_RFB, axis=1)

        data['EqSel_diff_from_opt'] = data.apply(self.get_diff_from_opt_EqSel, axis=1)
        
        return data


    def get_diff_from_opt_RFB(self, data):
        return data['max_grid_EL'] - max(data['RFB_maxContEL'], data['best_value_RFB'])
    
    def get_diff_from_opt_EqSel(self, data):
        return data['max_grid_EL'] - max(data['EqSel_maxContEL'], data['best_value_EqSel'])
    



    def _tidy_df(self, df):

        # avoid the "-1" case where has never failed:
        if df[df["fullDoseEL"]<=0].shape[0]:
            n_never_fail = df[df["fullDoseEL"]<=0].shape[0]
            print(f"{n_never_fail} runs never failed - need longer n_years. Filtering them out for now.")
            df = df[df["fullDoseEL"]>0]
        
        return df





    def _save_df(self):
        df_out = self.df
        
        filename = f"{self.folder}/par_scan/max_along_contour_{len(df_out)}.csv"

        print(f"Saving maximum along contour csv to: \n{filename}")

        df_out.to_csv(filename)


