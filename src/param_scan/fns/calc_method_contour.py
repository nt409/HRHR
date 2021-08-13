import pandas as pd
import numpy as np

from model.simulator import RunSingleTactic
from model.strategy_arrays import EqualResFreqBreakdownArray, \
    EqualSelectionArray


# TOC
# RunAlongContourDF

# ThisStratDetailedDF
# ThisStratSummaryDF




class RunAlongContourDFs:
    def __init__(self, rand_pars, grid_output,
                        n_cont_points, strat_name) -> None:
        
        print(f"Running contour method: {strat_name}")

        if strat_name=="RFB":
            strat_obj = EqualResFreqBreakdownArray(grid_output)
        elif strat_name=="EqSel":
            strat_obj = EqualSelectionArray(grid_output)
        else:
            raise Exception(f"invalid strat_name: {strat_name}")

        
        
        max_grid_EL = np.amax(grid_output.FY)

        self.df = ThisStratDetailedDF(rand_pars, 
                                    n_cont_points, 
                                    strat_obj,
                                    max_grid_EL).df
        
        
        self.summary = ThisStratSummaryDF(self.df, 
                                        strat_obj.name,
                                        strat_obj.level,
                                        max_grid_EL).df




# End of RunAlongContourDF






class ThisStratDetailedDF:
    def __init__(self, rand_pars, n_cont_points, strat_obj, max_grid_EL) -> None:
        
        self.rand_pars = rand_pars
        self.n_cont_points = n_cont_points

        self.strat_name = strat_obj.name
        self.level = strat_obj.level

        self.strat_obj = strat_obj
        self.max_grid_EL = max_grid_EL

        self.df = self._get_df()
        


    def _get_df(self):
        
        cntr = self.strat_obj.find_contours(self.rand_pars,
                                            self.n_cont_points)
        out = self._find_df(cntr)
        

        return out






    def _find_df(self, contour):
        
        df_this_run = pd.DataFrame()

        if not contour:
            return df_this_run


        for dose1, dose2, cont_quant in zip(contour['x'], contour['y'], contour['cont_vals']): 

            EL = self._get_EL_this_dose_combo(dose1, dose2)
            
            data = {"cont_quant": cont_quant,
                     "dose1": dose1,
                     "dose2": dose2,
                     "DS": dose1 + dose2,
                     "EL": EL}

            df_this_run = df_this_run.append(data, ignore_index=True)
        
        
        df_this_run['worked'] = max(df_this_run['EL'])>=self.max_grid_EL
        df_this_run['max_grid_EL'] = self.max_grid_EL

        return df_this_run





    def _get_EL_this_dose_combo(self, dose1, dose2):

        rp = self.rand_pars

        this_dose_conf = rp.get_single_conf(dose1, dose2)
                
        sing_run = RunSingleTactic(rp.fung_parms).run(this_dose_conf)
        
        EL = sing_run.failure_year
        
        return EL




# End of ThisStratDetailedDF












class ThisStratSummaryDF:

    def __init__(self, df, strat_name, level, max_grid_EL) -> None:
        self.df = df
        self.strat_name = strat_name
        self.level = level
        self.max_grid_EL = max_grid_EL

        self.LOW_THRESH = 1/3
        self.HIGH_THRESH = 2/3

        self.df = self.get_summary_df()





    def get_summary_df(self):
        df = self.df
        strat_name = self.strat_name

        if not df.shape[0]:
            return pd.DataFrame()

        lowDS_val = self._get_low_dose_max_EL()
        medDS_val = self._get_med_dose_max_EL()
        highDS_val = self._get_high_dose_max_EL()

        
        min_opt_dist_from_cntr = self._get_min_opt_dist_from_contour()
        

        data= {
            f"c_{strat_name[0]}_minContEL": min(df['EL']),
            f"c_{strat_name[0]}_maxContEL": max(df['EL']),

            f"c_{strat_name[0]}_min_opt_dist_from_contour": min_opt_dist_from_cntr,

            f"c_{strat_name[0]}_minDS": min(df['DS']),
            f"c_{strat_name[0]}_maxDS": max(df['DS']),
            
            f"c_{strat_name[0]}_worked_geq": self.max_grid_EL<=max(df['EL']),
            f"c_{strat_name[0]}_worked_equal": self.max_grid_EL==max(df['EL']),
            
            f"c_{strat_name[0]}_lowDoseMaxEL": lowDS_val,
            f"c_{strat_name[0]}_medDoseMaxEL": medDS_val,
            f"c_{strat_name[0]}_highDoseMaxEL": highDS_val,
            }

        return pd.DataFrame([data])


    


    
    
    def _get_low_dose_max_EL(self):
        df = self.df

        DS_thres_low = min(df['DS']) + self.LOW_THRESH*(max(df['DS']) - min(df['DS']))
        
        filt = df[df['DS'] < DS_thres_low]

        return self._get_max_if_df_non_empty(filt)
    





    def _get_med_dose_max_EL(self):
        df = self.df

        DS_thres_low = min(df['DS']) + self.LOW_THRESH*(max(df['DS']) - min(df['DS']))
        DS_thres_high = min(df['DS']) + self.HIGH_THRESH*(max(df['DS']) - min(df['DS']))

        filt = df[((df['DS'] >= DS_thres_low) & (df['DS'] <= DS_thres_high))]

        return self._get_max_if_df_non_empty(filt)

    




    def _get_high_dose_max_EL(self):
        df = self.df

        DS_thres_high = min(df['DS']) + self.HIGH_THRESH*(max(df['DS']) - min(df['DS']))

        filt = df[df['DS'] > DS_thres_high]
        
        return self._get_max_if_df_non_empty(filt)
        



    def _get_max_if_df_non_empty(self, df):
        if df.shape[0]:
            return max(df['EL'])
        else:
            return "NA"




    def _get_min_opt_dist_from_contour(self):
        df = self.df

        opt_df = df[df['EL']==self.max_grid_EL]

        vec = abs(opt_df["cont_quant"] - self.level)

        if len(vec):
            out = min(vec)
            return out
        else:
            return "NA"
        




# End of ThisStratSummaryDF





    




