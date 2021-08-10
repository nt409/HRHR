import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize



from utils.functions import RunSingleTactic, \
    logit10_difference, EqualResFreqBreakdownArray, \
    EqualSelectionArray

# TOC
# RunAlongContourDF
# DoseSumExtremesGetter
# ContourDoseFinder

# ThisStratSummaryDF

# find_FY_sel
# find_RFB


class RunAlongContourDFs:
    def __init__(self, rand_pars, grid_output,
                        n_points, strat_name) -> None:
        
        if strat_name=="RFB":
            level = 0
            strat_obj = EqualResFreqBreakdownArray(grid_output)
            
        elif strat_name=="EqSel":
            level = 0.5
            strat_obj = EqualSelectionArray(grid_output)
            
        else:
            raise Exception(f"invalid strat_name: {strat_name}")



        self.df = ThisStratDetailedDF(rand_pars, 
                                    n_points, 
                                    strat_name, 
                                    level,
                                    strat_obj).df
        
        
        max_grid_EL = np.amax(grid_output['FY'])

        self.summary = ThisStratSummaryDF(self.df, 
                                        strat_name,
                                        level,
                                        max_grid_EL).df










class ThisStratDetailedDF:
    def __init__(self, rand_pars, n_points, strat_name, 
                level, strat_obj) -> None:
        
        self.rand_pars = rand_pars
        self.n_points = n_points
        
        self.strat_name = strat_name
        self.level = level
        self.strat_obj = strat_obj

        self.df = self._get_df()


    def _get_df(self):

        cntr = self._find_cntr()        
        out = self._find_df(cntr)

        return out




    
    def _find_cntr(self):
        
        if not self.strat_obj.is_valid:
            print("strategy not valid:", self.strat_name)
            cntr_out = {}
            return cntr_out
        else:    
            DS_min_max = DoseSumExtremesGetter(self.strat_obj.array, self.level).min_max_DS
            
            cntr_out = ContourDoseFinder(self.rand_pars, self.strat_name,
                            DS_min_max, self.n_points, self.level).doses_out
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

        self.LOW_THRESH = 0.2
        self.HIGH_THRESH = 0.8

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





    






class DoseSumExtremesGetter:

    def __init__(self, z, level) -> None:

        cs = self._get_contour_set(z, level)
        
        self.min_max_DS = self._get_min_max_dose_sums(cs)


      


    def _get_contour_set(self, z, level):

        x, y = np.mgrid[0:1:z.shape[0]*1j, 0:1:z.shape[1]*1j]

        cs = plt.contour(x, y, z, levels=[level])
        
        return cs





    def _get_min_max_dose_sums(self, cs):

        if not cs.allsegs:
            print("Warning: cs.allsegs was empty!")
            # why not?
            # suspect too many 'nan's to create a proper contour
            return {}

        else:
            cont = cs.allsegs[0][0]
            
            x_vals = np.asarray(cont[:,0])
            y_vals = np.asarray(cont[:,1])

            dose_sum = x_vals + y_vals

            return dict(min=min(dose_sum), max=max(dose_sum))











class ContourDoseFinder:
    def __init__(self, rand_pars, strat_name, DS_min_max, n_points, level, tol=0.001) -> None:
        self.rand_pars = rand_pars
        
        if strat_name=="RFB":
            self.cont_quant_finder = find_RFB
            
        elif strat_name=="EqSel":
            self.cont_quant_finder = find_FY_sel
            
        else:
            raise Exception(f"invalid strat_name: {strat_name}")

        self.DS_min_max = DS_min_max
        self.n_points = n_points
        self.level = level

        self.CONT_DIST_THRESH = tol

        self.doses_out = self.get_doses_on_contour()
    



    def get_doses_on_contour(self):
        
        DS_bds = self.DS_min_max

        dose_sums = np.linspace(DS_bds['min'], DS_bds['max'], self.n_points)

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
                    print(self.DS_min_max)

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
    
    start_freqs = single_run['start_of_season']
    
    sr1 = start_freqs['RS'][1]/start_freqs['RS'][0]
    sr2 = start_freqs['SR'][1]/start_freqs['SR'][0]
    
    try:
        return sr1/(sr1+sr2)
    except Exception as e:
        print(f"FY finder warning/error: {e}, srs={(sr1,sr2)}")
        return None







def find_RFB(sing_run):
    
    fy_in = sing_run['failure_year']

    fy = int(fy_in)

    res_vecs = sing_run['res_vec_dict']
    
    rf1 = res_vecs['f1'][fy]
    
    rf2 = res_vecs['f2'][fy]

    try:
        return logit10_difference(rf1, rf2)
    except:
        print(f"problem with calculating delta_RFB for: {rf1, rf2}")
        return "NA"
