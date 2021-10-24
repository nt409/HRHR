import itertools
import numpy as np
from scipy.optimize import minimize
from math import floor

from model.utils import logit10_difference
from model.simulator import RunSingleTactic


# RFB Finder

def find_RFB_from_single_run(single_run):
    
    fy = single_run.failure_year

    end_freqs = single_run.end_freqs

    ef1 = end_freqs['RS']
    ef2 = end_freqs['SR']

    out = get_RFB(ef1, ef2, fy)
    return out




def get_RFB(end_freqs_1, end_freqs_2, fy):
    fy = int(fy)

    if not fy>1:
        return None
    
    r1 = end_freqs_1[fy-1]
    r2 = end_freqs_2[fy-1]

    try:
        return logit10_difference(r1, r2)
    except Exception as e:
        print(f"find_RFB warning/error: {e} \n rfs={r1, r2}")
        return None





# FY Sel Finder

def find_SFY_from_single_run(single_run):
    fy = single_run.failure_year
    
    start_freqs = single_run.start_freqs
    end_freqs = single_run.end_freqs

    ef1 = end_freqs['RS'][0]
    sf1 = start_freqs['RS'][0]
    
    ef2 = end_freqs['SR'][0]
    sf2 = start_freqs['SR'][0]

    out = get_FY_sel(ef1, sf1, ef2, sf2, fy)
    return out
    
    


def get_FY_sel(ef1, sf1, ef2, sf2, fy):
    """
    NB have changed so it compares end of season and start of season,
    not start of consecutive seasons. 
    
    This is because sexual reproduction can cause a change that wasn't due 
    to the tactic but was just down to ratios starting away from those expected
    in a perfectly mixed population.
    """
    
    fy = int(fy)

    if not fy>1:
        return None
    
    sr1 = ef1/sf1
    sr2 = ef2/sf2
    
    try:
        return sr1/(sr1+sr2)
    except Exception as e:
        print(f"find_FY_sel warning/error: {e} \n srs={(sr1,sr2)}")
        return None











class StrategyArray:
    def __init__(self, grid_output) -> None:
        self.grid_output = grid_output        
        self.contour_finder = ContourDoseFinder



    def _check_valid(self):
        """
        Check if Strategy is a possible tactic:
            - does the contour exist for this array?
        """
        return (np.nanmax(self.array)>self.level
                            and np.nanmin(self.array)<self.level)


    def find_contours(self, rand_pars, n_cont_points):

        if (not self.is_valid):
            return {}

        DS_extremes = DoseSumExtremes(self.array, self.level)

        if ((DS_extremes.min is None) or 
                (DS_extremes.max is None)):
            return {}
        

        cf = self.contour_finder(rand_pars,
                            self.name,
                            DS_extremes, 
                            n_cont_points,
                            self.level)
        
        contours = cf.get_doses_on_contour()
        
        return contours






class EqualResFreqBreakdownArray(StrategyArray):
    def __init__(self, grid_output) -> None:
        super().__init__(grid_output)

        self.level = 0
        self.name = "RFB"

        self.array = self._generate_array()
        self.is_valid = self._check_valid()
        
    
        
    def _generate_array(self):

        FYs = self.grid_output.FY
        end_freqs = self.grid_output.end_freqs_DA

        out = np.ones(FYs.shape)
        
        for i, j in itertools.product(range(out.shape[0]), 
                                        range(out.shape[1])):
            
            fy = FYs[i,j]
                
            ef1 = end_freqs['RS'][i,j,:]
            ef2 = end_freqs['SR'][i,j,:]

            out[i,j] = get_RFB(ef1, ef2, fy)

        return out




# * End of ERFB cls









class EqualSelectionArray(StrategyArray):
    def __init__(self, grid_output) -> None:
        super().__init__(grid_output)

        self.level = 0.5
        self.name = "EqSel"        

        self.array = self._generate_array()
        self.is_valid = self._check_valid()



    def _generate_array(self):

        FYs = self.grid_output.FY
        start_freqs = self.grid_output.start_freqs_DA
        end_freqs = self.grid_output.end_freqs_DA
        
        out = np.ones(FYs.shape)
        
        for i, j in itertools.product(range(out.shape[0]), 
                                        range(out.shape[1])):
            
            fy = FYs[i,j]

            ef1 = end_freqs['RS'][i,j,0]
            sf1 = start_freqs['RS'][i,j,0]
            
            ef2 = end_freqs['SR'][i,j,0]
            sf2 = start_freqs['SR'][i,j,0]

            out[i,j] = get_FY_sel(ef1, sf1, ef2, sf2, fy)

        return out





# * End of ES cls











class ContourDoseFinder:
    def __init__(self, rand_pars, 
                        strat_name,
                        DS_extremes,
                        n_cont_points,
                        level,
                        tol=0.001) -> None:

        self.rand_pars = rand_pars
        
        if strat_name=="RFB":
            self.cont_quant_finder = find_RFB_from_single_run
            
        elif strat_name=="EqSel":
            self.cont_quant_finder = find_SFY_from_single_run
            
        else:
            raise Exception(f"invalid strat_name: {strat_name}")

        self.DS_extremes = DS_extremes
        self.n_cont_points = n_cont_points
        self.level = level

        self.CONT_DIST_THRESH = tol



    def get_doses_on_contour(self):
        
        DS_bds = self.DS_extremes

        dose_sums = np.linspace(DS_bds.min, DS_bds.max, self.n_cont_points)

        x_list = []
        y_list = []
        cont_vals = []

        for ds in dose_sums:
            self.dose_sum = ds

            doses = self._get_doses_on_cntr_this_DS(ds)

            if doses is None:
                print("doses is None - means opt wasn't close or fy=1 - skipping to next ds")
                continue

            dist_from_contour = abs(self.model_cont_quant-self.level)
        
            if dist_from_contour < self.CONT_DIST_THRESH:
                x_list.append(doses['x'])
                y_list.append(doses['y'])
                cont_vals.append(self.model_cont_quant)
            else:
                print("\n")
                print("this run didn't get close to the contour?? ...")
                print("contour level:", self.model_cont_quant)
                print("dose sum:", self.dose_sum)
                print(f"{dist_from_contour=}")
                print(vars(DS_bds))


        return dict(x=x_list, y=y_list, cont_vals=cont_vals)




    
    def _get_doses_on_cntr_this_DS(self, ds):
        
        lower, upper = self._get_k_limits(ds)
        
        k0, bnds = self._get_k0_and_bds(lower, upper)

        if k0 is None:
            return None

        doses = self._get_doses(k0, bnds, ds)
        
        return doses






    def _get_k0_and_bds(self, lower, upper):
        """
        See if there is a point along the dose sum contour which gives FY>1
        - if so, return as k0
        - otherwise return None
        """

        ps = np.linspace(0,1,21)
        ks = []
        FYs = []

        for p in ps:
            k = (1-p)*lower + p*upper
            k_val, FY = self._valid_doses_this_value_k(k)
            
            ks.append(k_val)
            FYs.append(FY)
        
        if max(FYs)==1:
            return None, None
        
        FYs = np.asarray(FYs)
        inds = np.where(FYs==np.nanmax(FYs))[0]
        inds = list(inds)

        # if inds has several values, pick one in the middle
        middle_ind = inds[floor(len(inds)/2)]
        k = ks[middle_ind]

        
        ind_min = max(min(inds)-1, 0)
        ind_max = min(max(inds)+1, len(ks)-1)

        # limit search to near optimal values from basic search over k values
        k_min = ks[ind_min]
        k_max = ks[ind_max]

        return [k], ((k_min, k_max),)

            


    def _valid_doses_this_value_k(self, k):
        
        dose1 = 0.5*self.dose_sum - k
        dose2 = 0.5*self.dose_sum + k


        rp = self.rand_pars

        this_dose_conf = rp.get_single_conf(dose1, dose2)
        this_dose_conf.load_saved = False

        sing_run = RunSingleTactic(rp.fung_parms).run(this_dose_conf)

        FY = sing_run.failure_year
        
        return k, FY


    def _get_k_limits(self, ds):
        d = 0.5*ds

        # pick lower/upper so that each dose in [0,1]
        lower = max(-d, d - 1)
        upper = min( d,  1 - d)

        return lower, upper



    def _get_doses(self, k0, bnds, ds):
        # try:
        
        thisFit = minimize(self.objective_fn, k0, bounds=bnds)
        
        if self.model_cont_quant is None:
            print("model_cont_quant is None - return None")
            return None

        k = thisFit.x[0]

        doses = dict(
            x = 0.5*ds - k,
            y = 0.5*ds + k,
            )
        
        print(f"{doses=}")

        return doses

        # except Exception as e:
        #     print(e)
        #     return None
        





    def objective_fn(self, param):
        """
        1/ Updates self.model_cont_quant with RFB/SFY value here
        2/ Returns distance from contour if self.model_cont_quant!=None
        """

        k = param[0]

        dose1 = 0.5*self.dose_sum - k
        dose2 = 0.5*self.dose_sum + k

        rp = self.rand_pars

        this_dose_conf = rp.get_single_conf(dose1, dose2)
        this_dose_conf.load_saved = False

        sing_run = RunSingleTactic(rp.fung_parms).run(this_dose_conf)
        
        self.model_cont_quant = self.cont_quant_finder(sing_run)

        dist = self._get_distance_from_contour(k)

        # print(f"{dist=}")

        return dist




    def _get_distance_from_contour(self, k):
        """
        If quantity not defined (e.g. this dose combo is in unacceptable region)
        then return a high distance and one which pushes 
        the optimiser back towards the middle of the x+y=DS line
        """
        if self.model_cont_quant is None:
            return abs(k) + 2
        else:
            return (self.model_cont_quant - self.level)**2




                














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






