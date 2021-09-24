from model.strategy_arrays import EqualResFreqBreakdownArray
from param_scan.fns.calc_method_IVT import CheckStrategyUsingIVT_DF
import pandas as pd
import copy
import numpy as np


from model.simulator import RunGrid
from param_scan.fns.pars import RandomPars
from plotting.paper_figs import DoseSpaceScenarioSingle
from param_scan.fns.recalc_contour import RunAlongContourDFsReCalc

# TOC
# PostProcess
# ProcessedDF







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
        
        filename = f"{self.folder}/par_scan/paper/outcome_{df.shape[0]}.csv"
        print(f"Saving outcome data to: \n{filename}")
        out.to_csv(filename, index=False)
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

        df = copy.copy(self.processed_df)

        res = df.loc[:,[
                        "run",
                        "c_R_lowDoseMaxEL",
                        "c_R_medDoseMaxEL",
                        "c_R_highDoseMaxEL",
                        "max_grid_EL",
                        ]]

        res["low_opt"] = res["c_R_lowDoseMaxEL"] >= res["max_grid_EL"]
        res["med_opt"] = res["c_R_medDoseMaxEL"] >= res["max_grid_EL"]
        res["high_opt"] = res["c_R_highDoseMaxEL"] >= res["max_grid_EL"]

        # print(res.loc[res["run"]==328])
        # print(res.loc[res["high_opt"]])
        # runs = res.loc[res["high_opt"]]["run"]
        # x = self.par_df[self.par_df["run"].isin(runs)]
        # x = self.par_df[self.par_df["run"]==328]
        # x.sort_values(by="sr_prop", inplace=True)
        # print(x)
        # exit()
        # exit()
        
        out = pd.DataFrame()

        for row in ["low", "med", "high"]:
            data = dict(count = sum(res[row+"_opt"]),
                        pc = 100*sum(res[row+"_opt"])/res.shape[0],
                        dose = row)

            out = out.append(data, ignore_index=True)

        print(out)
        
        filename = f"{self.folder}/par_scan/paper/high_or_low_dose_{res.shape[0]}.csv"
        print(f"Saving high or low dose csv to: \n{filename}")
        out.to_csv(filename)




  
    
    
    def re_run_grid(self, NDoses, run_indices, plot=True):

        df_test = self.get_params_for_specific_runs(run_indices)

        for ii in range(df_test.shape[0]):

            pars = df_test.iloc[int(ii),:]
           
            this_run_ind = int(df_test.iloc[int(ii),:].run)

            print(f"\nRe-running run: {this_run_ind} \n")

            rp = self._get_RPs(pars, NDoses)
            
            grid_output = RunGrid(rp.fung_parms).run(rp.grid_conf)

            conf_str = rp.grid_conf.config_string_img

            FY = grid_output.FY
            opt_region = FY == np.amax(FY)
                
            n_opt_doses = opt_region.sum()

            print(f"Number of optimal dose combos: {n_opt_doses}")

            if plot:
                conf_str = conf_str.replace("param_scan/", f"param_scan/run={this_run_ind}_")
                DoseSpaceScenarioSingle(grid_output, conf_str)
        
        # return last one calculated
        return grid_output




    def re_run_cont(self, NDoses, N_cont_doses, DS_lim, run_indices):

        df_test = self.get_params_for_specific_runs(run_indices)

        out = pd.DataFrame()

        for ii in range(df_test.shape[0]):

            pars = df_test.iloc[int(ii),:]

            this_run_ind = int(df_test.iloc[int(ii),:].run)

            print(f"\nRe-running run: {this_run_ind} \n")

            rp = self._get_RPs(pars, NDoses)

            grid_output = RunGrid(rp.fung_parms).run(rp.grid_conf)

            # print(np.where(np.array(grid_output.FY)==np.amax(grid_output.FY)))

            RFB_dfs = RunAlongContourDFsReCalc(rp, grid_output, N_cont_doses, DS_lim, "RFB")

            data = dict(run=this_run_ind,
                    c_R_maxContEL=RFB_dfs.summary.c_R_maxContEL,
                    max_grid_EL=np.amax(grid_output.FY))
            
            out = out.append(data, ignore_index=True)

            # plot output
            conf_str = rp.grid_conf.config_string_img
            conf_str = conf_str.replace("param_scan/", f"param_scan/run={this_run_ind}_")

            DoseSpaceScenarioSingle(grid_output, conf_str)

        for col in ["c_R_maxContEL", "max_grid_EL", "run"]:
            out[col] = out[col].astype("int")
        out["worked"] = out["c_R_maxContEL"] >= out["max_grid_EL"]


        ind_str = ",".join([str(rr) for rr in run_indices])
        filename = f"{self.folder}/par_scan/re_run/cont_{NDoses}_{N_cont_doses}_{ind_str}.csv"
        print(f"Saving re-run to: {filename}")
        print(out)
        out.to_csv(filename, index=False)



    def re_run_IVT(self, NDoses, run_indices):

        df_test = self.get_params_for_specific_runs(run_indices)

        out = pd.DataFrame()

        for ii in range(df_test.shape[0]):

            pars = df_test.iloc[int(ii),:]

            this_run_ind = int(df_test.iloc[int(ii),:].run)

            print(f"\nRe-running run: {this_run_ind} \n")

            rp = self._get_RPs(pars, NDoses)

            grid_output = RunGrid(rp.fung_parms).run(rp.grid_conf)

            ##
            # print(pars)
            # exit()

            mask = np.where(grid_output.FY==12)



            # for ii in range(len(mask[0])):
            #     print("\n", ii)
            #     yy = grid_output.yield_array[mask[0][ii], mask[1][ii], 11]
            #     print(yy)

            #     # print(grid_output.start_freqs_DA["SR"][mask[0][ii], mask[1][ii],12] > grid_output.start_freqs_DA["RS"][mask[0][ii], mask[1][ii],12])
            #     print(grid_output.start_freqs_DA["SR"][mask[0][ii], mask[1][ii],12])
            #     print(grid_output.start_freqs_DA["RS"][mask[0][ii], mask[1][ii],12])
            #     print(grid_output.start_freqs_DA["RR"][mask[0][ii], mask[1][ii],12])
            #     # print(grid_output.end_freqs_DA["SR"][mask[0][ii], mask[1][ii],12] > grid_output.start_freqs_DA["RS"][mask[0][ii], mask[1][ii],12])

            
            # z = EqualResFreqBreakdownArray(grid_output).array
            # print(z[mask])
            # ##
            # exit()



            x = CheckStrategyUsingIVT_DF(grid_output, "RFB")

            data = dict(run=this_run_ind,
                    IVT_bv=x.best_value,
                    max_grid_EL=np.amax(grid_output.FY))

            out = out.append(data, ignore_index=True)

            # plot output
            conf_str = rp.grid_conf.config_string_img
            conf_str = conf_str.replace("param_scan/", f"param_scan/run={this_run_ind}_")

            DoseSpaceScenarioSingle(grid_output, conf_str)

        for col in ["IVT_bv", "max_grid_EL", "run"]:
            out[col] = out[col].astype("int")
        out["worked"] = out["IVT_bv"] >= out["max_grid_EL"]

        ind_str = ",".join([str(rr) for rr in run_indices])
        filename = f"{self.folder}/par_scan/re_run/IVT_{NDoses}_{ind_str}.csv"
        print(f"Saving re-run to: {filename}")
        print(out)
        out.to_csv(filename, index=False)





    @staticmethod
    def _get_RPs(pars, NDoses):

        config = {'load_saved': True, 'save': True, 'n_years': 35}

        RP = RandomPars(config, None)

        RP.get_inoc_dict(pars["RS"], pars["SR"], pars["RR"])
        
        RP.fung_parms = RP.get_fung_parms_dict(pars["omega_1"], pars["omega_2"], 
                                pars["delta_1"], pars["delta_2"], 
                                pars["theta_1"], pars["theta_2"])
        
        RP.sr_prop = pars["sr_prop"]

        RP.path_and_fung_pars = (pars["RS"], pars["SR"], pars["RR"], 
                    pars["omega_1"], pars["omega_2"], 
                    pars["delta_1"], pars["delta_2"], 
                    pars["theta_1"], pars["theta_2"])

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

        # data['ERFB_diff_from_opt'] = data.apply(self.get_diff_from_opt_RFB, axis=1)
        # data['ESFY_diff_from_opt'] = data.apply(self.get_diff_from_opt_EqSel, axis=1)

        data['ERFB_diff_from_opt'] = data["c_R_maxContEL"] - data["max_grid_EL"]
        data['ESFY_diff_from_opt'] = data["c_E_maxContEL"] - data["max_grid_EL"]
        
        data['ESFYL_diff_from_opt'] = data["c_E_lowDoseMaxEL"] - data["max_grid_EL"]


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



