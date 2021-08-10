import pandas as pd
from scipy.sparse.construct import rand

from utils.functions import RunGrid

from .pars import RandomPars
from .calc_method_contour import RunAlongContourDFs
from .calc_method_IVT import CheckStrategyUsingIVT_DF
from .calc_other_strats import OtherStratsDF






class SinglePSRun:
    def __init__(self, config, run_index) -> None:
        
        self.config = config
        
        self.run_index = run_index

        rand_pars = self._get_rand_pars_obj()
        
        grid_output = RunGrid(rand_pars.fung_parms).run(rand_pars.grid_conf)
        
        self.output = self._get_output(rand_pars, grid_output)



    def _get_output(self, rand_pars, grid_output):
        
        par_df = rand_pars.par_df
        
        RFB_dfs = RunAlongContourDFs(rand_pars, grid_output,
                        self.config['n_cont_points'], "RFB")

        EqSel_dfs = RunAlongContourDFs(rand_pars, grid_output,
                        self.config['n_cont_points'], "EqSel")

        other_strats_df = OtherStratsDF(grid_output).df

        IVT_RFB_df = CheckStrategyUsingIVT_DF(grid_output, "RFB").df

        IVT_EqSel_df = CheckStrategyUsingIVT_DF(grid_output, "EqSel").df

        summary_df = pd.concat([ 
                    RFB_dfs.summary, 
                    EqSel_dfs.summary,
                    other_strats_df,
                    IVT_RFB_df,
                    IVT_EqSel_df],
                    axis=1)


        summary_df['run'] = self.run_index
        RFB_dfs.df['run'] = self.run_index
        EqSel_dfs.df['run'] = self.run_index

        out = ScanOutput()
        out.add_data(par_df, summary_df, RFB_dfs.df, EqSel_dfs.df)

        return out
        



    def _get_rand_pars_obj(self):
    
        RP = RandomPars(self.config, self.run_index)

        RP.find_pars()
        
        RP.get_all_parms_dict()

        RP.get_grid_conf(self.config["grid_number"])

        return RP









class ScanOutput:
    """
    Output for one or many PS runs.
    """

    def __init__(self) -> None:
        self.par_df = pd.DataFrame()
        self.summary_df = pd.DataFrame()
        self.RFB_df = pd.DataFrame()
        self.EqSel_df = pd.DataFrame()


    def add_data(self, par_df, summary_df, RFB_df, EqSel_df):
        self.par_df = par_df
        self.summary_df = summary_df
        self.RFB_df = RFB_df
        self.EqSel_df = EqSel_df


    def add_new_output(self, scan_output):
        self.par_df = pd.concat([self.par_df, scan_output.par_df], axis=0)
        self.summary_df = pd.concat([self.summary_df, scan_output.summary_df], axis=0)
        self.RFB_df = pd.concat([self.RFB_df, scan_output.RFB_df], axis=0)
        self.EqSel_df = pd.concat([self.EqSel_df, scan_output.EqSel_df], axis=0)


    def save(self, config, seed):
        
        par_str = config['par_str']

        for key, df in vars(self).items():

            filename = (config['folder_save'] + "/par_scan/" +
                                f"{key}_seed={seed}_{par_str}.csv")
            
            print(f"\n Random Scan {key}, saved as:\n {filename}")

            df_out = self._move_run_to_first_col(df)

            df_out.to_csv(filename, index=False)




    def _move_run_to_first_col(self, df):
        cols = list(df.columns)
        new_cols = ['run'] + [cc for cc in cols if cc!='run']
        out = df.loc[:, new_cols]
        return out
