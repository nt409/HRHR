import pandas as pd
import numpy as np
from scipy.optimize import minimize

from model.simulator import RunSingleTactic
from model.utils import logit10_difference, \
    EqualResFreqBreakdownArray, \
    EqualSelectionArray


# TOC
# RunAlongContourDF
# DoseSumExtremes
# ContourDoseFinder

# ThisStratSummaryDF

# find_FY_sel
# find_RFB


class RunAlongContourDFs:
    def __init__(self, rand_pars, grid_output,
                        n_cont_points, strat_name) -> None:
        
        print(f"Running contour method: {strat_name}")

        if strat_name=="RFB":
            level = 0
            strat_obj = EqualResFreqBreakdownArray(grid_output)
            
        elif strat_name=="EqSel":
            level = 0.5
            strat_obj = EqualSelectionArray(grid_output)
            
        else:
            raise Exception(f"invalid strat_name: {strat_name}")

        
        
        max_grid_EL = np.amax(grid_output['FY'])

        self.df = ThisStratDetailedDF(rand_pars, 
                                    n_cont_points, 
                                    strat_name, 
                                    level,
                                    strat_obj,
                                    max_grid_EL).df
        
        
        self.summary = ThisStratSummaryDF(self.df, 
                                        strat_name,
                                        level,
                                        max_grid_EL).df










class ThisStratDetailedDF:
    def __init__(self, rand_pars, n_cont_points, strat_name, 
                level, strat_obj, max_grid_EL) -> None:
        
        self.rand_pars = rand_pars
        self.n_cont_points = n_cont_points
        self.strat_name = strat_name
        self.level = level
        self.strat_obj = strat_obj
        self.max_grid_EL = max_grid_EL

        self.df = self._get_df()


    def _get_df(self):

        cntr = self._find_cntr()        
        out = self._find_df(cntr)
        out['worked'] = max(out['EL'])>=self.max_grid_EL
        out['max_grid_EL'] = self.max_grid_EL

        return out




    
    def _find_cntr(self):
        
        if not self.strat_obj.is_valid:
            print("strategy not valid:", self.strat_name)
            return {}
        else:    
            DS_extremes = DoseSumExtremes(self.strat_obj.array, self.level)
            
            if DS_extremes.min is None or DS_extremes.max is None:
                return {}
            
            cntr_out = ContourDoseFinder(self.rand_pars, self.strat_name,
                            DS_extremes, self.n_cont_points, self.level).doses_out
            return cntr_out






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
        
        
        return df_this_run





    def _get_EL_this_dose_combo(self, dose1, dose2):

        rp = self.rand_pars

        this_dose_conf = rp.get_single_conf(dose1, dose2)
                
        sing_run = RunSingleTactic(rp.fung_parms).run(this_dose_conf)
        
        EL = sing_run['failure_year']
        
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





    






class DoseSumExtremes:

    def __init__(self, z, level) -> None:
        self.min = None
        self.max = None
        
        self.get_min_and_max_valid_DS(z, level)



    def get_min_and_max_valid_DS(self, z, level):
        """
        Check for each dose sum whether there are values above and below 'level'
        """
        # relative to 0
        zz = np.array(z) - level

        ds_vec = np.linspace(0, 2, -1+2*z.shape[0])

        valid_ds_list = []

        for ds_ind in range(len(ds_vec)):
            is_valid = self._check_this_ds_straddles_level(zz, ds_ind, z.shape[0])
            
            if is_valid:
                valid_ds_list.append(ds_vec[ds_ind])


        if not len(valid_ds_list):
            return None

        self.min = min(valid_ds_list)
        self.max = max(valid_ds_list)






    def _check_this_ds_straddles_level(self, zz, ds_ind, n):
        vals = []

        bottom = max(0, 1 + ds_ind-n)
        top = min(ds_ind, n-1)
        
        for ii in range(bottom, 1+top):
            if not np.isnan(zz[ii, ds_ind-ii]):
                vals.append(zz[ii, ds_ind-ii])

        if not len(vals):
            return False

        if min(vals)<0 and max(vals)>0:
            return True
        
        return False












class ContourDoseFinder:
    def __init__(self, rand_pars, strat_name, DS_extremes, n_cont_points, level, tol=0.001) -> None:
        self.rand_pars = rand_pars
        
        if strat_name=="RFB":
            self.cont_quant_finder = find_RFB
            
        elif strat_name=="EqSel":
            self.cont_quant_finder = find_FY_sel
            
        else:
            raise Exception(f"invalid strat_name: {strat_name}")

        self.DS_extremes = DS_extremes
        self.n_cont_points = n_cont_points
        self.level = level

        self.CONT_DIST_THRESH = tol

        self.doses_out = self.get_doses_on_contour()
    



    def get_doses_on_contour(self):
        
        DS_bds = self.DS_extremes

        dose_sums = np.linspace(DS_bds.min, DS_bds.max, self.n_cont_points)

        x_list = []
        y_list = []
        cont_vals = []

        for ds in dose_sums:
            self.dose_sum = ds

            k = self._get_doses_on_contour_single_DS(ds)

            if k is None:
                continue

            
            try:

                # only add if have got close to the contour (but why doesn't it get close otherwise?)
                
                # if self.model_cont_quant is NA then this doesn't work

                dist_from_contour = abs(self.model_cont_quant-self.level)
            
                if dist_from_contour < self.CONT_DIST_THRESH:
                    x_list.append(0.5*ds - k)
                    y_list.append(0.5*ds + k)
                    cont_vals.append(self.model_cont_quant)
                else:
                    print("\n")
                    print("this run didn't get close to the contour?? ...")
                    print("contour level:", self.model_cont_quant)
                    print("dose sum:", self.dose_sum)
                    print(DS_bds)

            except Exception as e:
                print(e)

                

        return dict(x=x_list, y=y_list, cont_vals=cont_vals)




    
    def _get_doses_on_contour_single_DS(self, ds):
        
        lower, upper = self._get_bnds(ds)
        
        x0 = [0.5*(lower+upper)]
        
        bnds = ((lower, upper), )

        try:
            thisFit = minimize(self.objective_fn, x0, bounds=bnds)
            return thisFit.x[0]
        except Exception as e:
            print(e)
            return None
        




    def _get_bnds(self, ds):
        d = 0.5*ds

        # pick lower/upper so that each dose in [0,1]
        lower = max(-d, d - 1)
        upper = min( d,  1 - d)

        return lower, upper






    def objective_fn(self, param):

        k = param[0]

        dose1 = 0.5*self.dose_sum - k
        dose2 = 0.5*self.dose_sum + k

        rp = self.rand_pars

        this_dose_conf = rp.get_single_conf(dose1, dose2)
        this_dose_conf.load_saved = False

        sing_run = RunSingleTactic(rp.fung_parms).run(this_dose_conf)
        
        self.model_cont_quant = self.cont_quant_finder(sing_run)

        dist = self._get_distance_from_contour(k)

        return dist




    def _get_distance_from_contour(self, k):
        """
        If quantity not defined (e.g. this dose combo is in unacceptable region)
        then return a high distance and one which pushes 
        the optimiser back towards the middle of the x+y=DS line
        """

        try:
            return (self.model_cont_quant - self.level)**2
        except:
            return abs(k) + 5




                












def find_FY_sel(single_run):
    
    start_freqs = single_run['start_freqs']
    end_freqs = single_run['end_freqs']
    
    sr1 = end_freqs['RS'][0]/start_freqs['RS'][0]
    sr2 = end_freqs['SR'][0]/start_freqs['SR'][0]
    
    try:
        return sr1/(sr1+sr2)
    except Exception as e:
        print(f"find_FY_sel warning/error: {e} \n srs={(sr1,sr2)}")
        return None







def find_RFB(sing_run):
    
    fy = sing_run['failure_year']
    fy = int(fy)

    end_freqs = sing_run['end_freqs']
    
    r1 = end_freqs['RS'][fy-1]
    r2 = end_freqs['SR'][fy-1]

    try:
        return logit10_difference(r1, r2)
    except Exception as e:
        print(f"find_RFB warning/error: {e} \n rfs={r1, r2}")
        return "NA"
