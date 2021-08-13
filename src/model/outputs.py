import numpy as np




class ODEStates:
    def __init__(self, t) -> None:

        self.t = t

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


    def update_y(self, y):

        ind = self.index

        length = y.shape[1]

        self.S[ind:ind+length] = y[0, :]
        
        self.ERR[ind:ind+length] = y[1, :]
        self.ERS[ind:ind+length] = y[2, :]
        self.ESR[ind:ind+length] = y[3, :]
        self.ESS[ind:ind+length] = y[4, :]
        
        self.IRR[ind:ind+length] = y[5, :]
        self.IRS[ind:ind+length] = y[6, :]
        self.ISR[ind:ind+length] = y[7, :]
        self.ISS[ind:ind+length] = y[8, :]
        
        self.R[ind:ind+length] = y[9, :]
        
        self.PRR[ind:ind+length] = y[10, :]
        self.PRS[ind:ind+length] = y[11, :]
        self.PSR[ind:ind+length] = y[12, :]
        self.PSS[ind:ind+length] = y[13, :]

        self.fung_1[ind:ind+length] = y[14, :]
        self.fung_2[ind:ind+length] = y[15, :]

        self.index = ind + length
    
    def delete_unnecessary_vars(self):
        delattr(self, "index")










class SimOutput:
    def __init__(self, t_vec) -> None:
        self.states = ODEStates(t_vec)

        self.final_res_vec_dict = None
        self.end_freqs = None
        self.selection = None
        self.yield_val = None











class SingleTacticOutput:
    def __init__(self, yield_thresh, res_props_in, strain_names, n_years, df_yield) -> None:

        self.yield_thresh = yield_thresh
        self.n_years = n_years
        self.df_yield = df_yield
        self.strain_names = strain_names
        

        self.failure_year = 0
        
        self.yield_vec = np.zeros(n_years)

        self.res_vec_dict = self._init_res_vec_dict(res_props_in)

        self.selection_vec_dict = self._init_dict_of_zeros_vecs(['f1', 'f2'], n_years+1)
        
        # post-sex from previous year, ready for start of season
        self.start_freqs = self._init_dict_of_zeros_vecs(strain_names, n_years+1)

        # end of season, pre-sex
        self.end_freqs = self._init_dict_of_zeros_vecs(strain_names, n_years+1)

        self.states_list = []



    def _init_dict_of_zeros_vecs(self, keys, length):
        out = {}
        
        for key in keys:
            out[key] = np.zeros(length)

        return out
    

    def _init_res_vec_dict(self, res_props):
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
        
        if ((self.yield_vec[yr]<self.yield_thresh) and 
                (self.failure_year==0)):
            self.failure_year = yr+1


    def delete_unnecessary_vars(self):
        delattr(self, "yield_thresh")
        # delattr(self, "conf")
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
