import numpy as np
from math import floor




class ODEStates:
    def __init__(self, params, tvs) -> None:

        self.params = params

        self.t = self._get_t(tvs)

        n_points = len(self.t)
        
        # start filling up from 0th index
        self.index = 0

        # state values
        self.S = np.zeros(n_points)
        
        self.ERR = np.zeros(n_points)
        self.ERS = np.zeros(n_points)
        self.ESR = np.zeros(n_points)
        self.ESS = np.zeros(n_points)

        self.IRR = np.zeros(n_points)
        self.IRS = np.zeros(n_points)
        self.ISR = np.zeros(n_points)
        self.ISS = np.zeros(n_points)

        self.R = np.zeros(n_points)

        self.PRR = np.zeros(n_points)
        self.PRS = np.zeros(n_points)
        self.PSR = np.zeros(n_points)
        self.PSS = np.zeros(n_points)

        self.fung_1 = np.zeros(n_points)
        self.fung_2 = np.zeros(n_points)


    
 
    @staticmethod    
    def _get_t(tvs):

        out = np.concatenate([tvs[ii][:-1] if ii!=3 else tvs[ii]
                                for ii in range(len(tvs))])
        return out



    def update_y(self, y):
        p = self.params

        ind = self.index

        length = y.shape[1]

        self.S[ind:ind+length] = y[p.S_ind, :]
        
        self.ERR[ind:ind+length] = y[p.ERR_ind, :]
        self.ERS[ind:ind+length] = y[p.ERS_ind, :]
        self.ESR[ind:ind+length] = y[p.ESR_ind, :]
        self.ESS[ind:ind+length] = y[p.ESS_ind, :]
        
        self.IRR[ind:ind+length] = y[p.IRR_ind, :]
        self.IRS[ind:ind+length] = y[p.IRS_ind, :]
        self.ISR[ind:ind+length] = y[p.ISR_ind, :]
        self.ISS[ind:ind+length] = y[p.ISS_ind, :]
        
        self.R[ind:ind+length] = y[p.R_ind, :]
        
        self.PRR[ind:ind+length] = y[p.PRR_ind, :]
        self.PRS[ind:ind+length] = y[p.PRS_ind, :]
        self.PSR[ind:ind+length] = y[p.PSR_ind, :]
        self.PSS[ind:ind+length] = y[p.PSS_ind, :]

        self.fung_1[ind:ind+length] = y[p.fung_1_ind, :]
        self.fung_2[ind:ind+length] = y[p.fung_2_ind, :]

        self.index = ind + length
    














class SimOutput:
    def __init__(self, params) -> None:

        self.params = params

        self.seg_times = [self.params.T_emerge,
                        self.params.T_GS32,
                        self.params.T_GS39,
                        self.params.T_GS61,
                        self.params.T_GS87]
        
        self.seg_names = ["start", "spray_1", "spray_2", "yield"]

        self.t_vecs = self._get_list_of_time_vecs()

        self.states = ODEStates(params, self.t_vecs)




    def _get_list_of_time_vecs(self):
        
        seg_ts = self.seg_times
        
        sum_ns = 0

        list_of_tvs = []

        for ii, segment in enumerate(self.seg_names):

            if segment=="yield":
                # makes sure total number of points is self.params.t_points
                n = 3 + self.params.t_points - sum_ns

            else:
                # make n so that values are approx self.params.dt apart
                n = 1 + (seg_ts[ii+1]-seg_ts[ii])/self.params.dt
                n = floor(n)
                sum_ns += n

            time_vec = np.linspace(seg_ts[ii], seg_ts[ii+1], n)

            if segment=="yield":
                self.t_yield = time_vec

            list_of_tvs.append(time_vec)

        return list_of_tvs




    def get_yield_contributing_y(self, y):
        p = self.params

        self.y_yield = (y[p.S_ind,:] + 
                            y[p.ERR_ind,:] +
                            y[p.ERS_ind,:] +
                            y[p.ESR_ind,:] +
                            y[p.ESS_ind,:])






    def delete_unnecessary_vars(self):
        delattr(self.states, "index")
        delattr(self.states, "params")
        
        delattr(self, "params")
        delattr(self, "seg_times")
        delattr(self, "seg_names")
        delattr(self, "t_vecs")
        delattr(self, "t_yield")
        delattr(self, "y_yield")
















class SingleTacticOutput:
    def __init__(self, params, config, strain_names, n_years, df_yield) -> None:

        self.params = params
        self.conf = config
        self.n_years = n_years
        self.df_yield = df_yield
        self.strain_names = strain_names
        

        self.failure_year = 0
        
        self.yield_vec = np.zeros(n_years)

        self.res_vec_dict = self._initialise_res_vec_dict(self.conf.res_props)

        self.selection_vec_dict = self._initialise_dict_of_vecs(['f1', 'f2'], n_years+1)
        
        # post-sex from previous year, ready for start of season
        self.start_freqs = self._initialise_dict_of_vecs(strain_names, n_years+1)

        # end of season, pre-sex
        self.end_freqs = self._initialise_dict_of_vecs(strain_names, n_years+1)

        self.states_list = []



    def _initialise_dict_of_vecs(self, keys, length):
        out = {}
        
        for key in keys:
            out[key] = np.zeros(length)

        return out
    

    def _initialise_res_vec_dict(self, res_props):
        out = {}
        keys = ['f1', 'f2']
        
        for key in keys:
            out[key] = np.zeros(self.n_years+1)
            # set first year
            out[key][0] = res_props[key]

        return out






    def add_new_sim_output(self, sim_out, yr):
        
        self.yield_vec[yr] = 100*(sim_out.yield_val/self.df_yield)

        self._update_end_freqs(sim_out, yr)

        self._update_selection_vec_dict(sim_out, yr+1)
        
        self._update_res_vec_dict(sim_out, yr+1)
        
        self._update_failure_year(yr)
        
        self.states_list.append(sim_out.states)






    
    def _update_selection_vec_dict(self, sim_out, yr):
        for key in ['f1', 'f2']:
            self.selection_vec_dict[key][yr] = sim_out.selection[key]


    def _update_res_vec_dict(self, sim_out, yr):
        for key in ['f1', 'f2']:
            self.res_vec_dict[key][yr] = sim_out.final_res_vec_dict[key]

    
    def update_start_freqs(self, values, yr):
        for key in self.strain_names:
            self.start_freqs[key][yr] = values[key]
    

    def _update_end_freqs(self, sim_out, yr):
        for key in self.end_freqs.keys():
            self.end_freqs[key][yr] = sim_out.end_freqs[key]


    def _update_failure_year(self, yr):
        """
        Set failure year if:
        - yield is below threshold
        - is first time it has dropped below threshold
        """
        
        if ((self.yield_vec[yr]<self.params.yield_threshold) and 
                (self.failure_year==0)):
            self.failure_year = yr+1


    def delete_unnecessary_vars(self):
        delattr(self, "params")
        delattr(self, "conf")
        delattr(self, "n_years")
        delattr(self, "df_yield")
        delattr(self, "strain_names")









class GridTacticOutput:
    def __init__(self, n_doses, n_years) -> None:

        self.LTY = np.zeros((n_doses, n_doses))
        self.TY = np.zeros((n_doses, n_doses))
        self.FY = np.zeros((n_doses, n_doses))
        
        self.yield_array = np.zeros((n_doses, n_doses, n_years))
        
        fung_keys = ['f1', 'f2']
        self.selection_DA = self._get_dict_of_zero_arrays(fung_keys, (n_doses, n_doses, n_years+1))
        self.res_vec_DA = self._get_dict_of_zero_arrays(fung_keys, (n_doses, n_doses, n_years+1))

        strain_keys = ['RR', 'RS', 'SR', 'SS']
        self.start_freqs_DA = self._get_dict_of_zero_arrays(strain_keys, (n_doses, n_doses, n_years+1))
        self.end_freqs_DA = self._get_dict_of_zero_arrays(strain_keys, (n_doses, n_doses, n_years+1))
    
    
    
    @staticmethod
    def _get_dict_of_zero_arrays(keys, shape):
        out = {}
        for key in keys:
            out[key] = np.zeros(shape)
        return out




    def update_dicts_of_arrays(self, data, f1_ind, f2_ind):
        self.selection_DA = self._update_dict_array_this_dose(
                                            self.selection_DA, data, 
                                            f1_ind, f2_ind, "selection_vec_dict")

        self.res_vec_DA = self._update_dict_array_this_dose(
                                            self.res_vec_DA, data, 
                                            f1_ind, f2_ind, "res_vec_dict")

        self.start_freqs_DA = self._update_dict_array_this_dose(
                                            self.start_freqs_DA, data, 
                                            f1_ind, f2_ind, "start_freqs")

        self.end_freqs_DA = self._update_dict_array_this_dose(
                                            self.end_freqs_DA, data, 
                                            f1_ind, f2_ind, "end_freqs")



    def _update_dict_array_this_dose(self, to_update, data, f1_ind, f2_ind, attr):
        
        for key_ in to_update.keys():
            calculated = vars(data)[attr]
            to_update[key_][f1_ind,f2_ind,:] = calculated[key_]
        
        return to_update
