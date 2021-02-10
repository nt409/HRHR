import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pickle
import copy
import json
from scipy.interpolate import RegularGridInterpolator
from math import ceil, log, floor, log10, exp
import warnings
##
# import sys
# sys.path.insert(0, '/Functions_and_plotting/')
##
from utils.functions_HRHR import master_loop_grid_of_tactics
from utils.Optimal_simulator_functions_and_animator import FD_space
from utils.plotter_HRHR    import Overlay_plotter, mosaic_plot_2d, cube_scatter_plot, surface_plot_3d, generate_array, dose_choice_plotter, phi_space_navigator
from utils.parameters_HRHR import params, params_dict
#----------------------------------------------------------------------------------------------
def recursion(F_array,f_no,n_d,dose_array_LTY,LTY_A2,LTY_int,log_phi_min = -5,cube_arg=0):
    n_p     = F_array.shape[0]
    phi_vec = np.linspace(log_phi_min,0,n_p)
    x,y,z   = [phi_vec for i in range(3)]
    Int_F   = RegularGridInterpolator((x,y,z),F_array,bounds_error = False,fill_value=f_no)
    ##
    prr2, prs2, psr2, Yield = [2*np.ones((n_p,n_p,n_p,n_d,n_d)) for ii in range(4)]
    LTY_A1      = np.zeros((n_p,n_p,n_p,n_d,n_d))
    LTY_A2      = copy.deepcopy(LTY_A2) # since interpolation function was changing when the array did.
    for i2 in range(n_p):
        for j2 in range(n_p):
            for k2 in range(n_p):
                i = (n_p-1)-i2 # counting backwards allows us to calculate minimal numbers assuming 'strictly better points' well behaved using cube argument below
                j = (n_p-1)-j2
                k = (n_p-1)-k2
                if F_array[i,j,k]==f_no:
                    prr = 10**(phi_vec[i])
                    prs = 10**(phi_vec[j])
                    psr = 10**(phi_vec[k])
                    pss = 1 - prr - prs - psr
                    if pss>0:
                        output = master_loop_grid_of_tactics(n_d,1,p_rr=prr,p_rs=prs,p_sr=psr,p_ss=pss)
                        prr2[i,j,k,:,:]  = output['PRR_array'][:,:,1] # only one season
                        prs2[i,j,k,:,:]  = output['PRS_array'][:,:,1] # only one season
                        psr2[i,j,k,:,:]  = output['PSR_array'][:,:,1] # only one season
                        Yield[i,j,k,:,:] = output['Yield'][:,:,0]     # only one season
                        ##
                        # update F_array then LTY_A1
                        for d1 in range(n_d):
                            for d2 in range(n_d):
                                X = [log10(prr2[i,j,k,d1,d2]),log10(prs2[i,j,k,d1,d2]),log10(psr2[i,j,k,d1,d2])]
                                if Yield[i,j,k,d1,d2]>params.Yield_threshold and F_array[i,j,k]==f_no: # so that haven't already been 'upgraded'
                                    # print(f_no,Int_F(X),X,(i,j,k,d1,d2))
                                    if min(X)<log_phi_min:
                                        warnings.warn('Warning, fill value used for X =',X) # print
                                    if Int_F(X) > f_no - 0.05:# and Int_F(X)< f_no+1: # if we get sent to somewhere with at least as good a f_no # same f_no
                                        F_array[i,j,k] = f_no + 1
                                        # print(F_array[i,:,:],(i,j,k))
                                        ### cube argument - aiming to reduce computational time. Probably no longer applicable
                                        if cube_arg !=0:
                                            for i3 in range(i):
                                                for j3 in range(j):
                                                    for k3 in range(k):
                                                        print(i3,j3,k3,'cube argument',i,j,k)
                                                        F_array[i3,j3,k3] = f_no + 1
                                        ####################
                                if Yield[i,j,k,d1,d2]>params.Yield_threshold and max(X)<0 and min(X)>log_phi_min and F_array[i,j,k]>=f_no: # !! was ==f_no # !! this excludes some potentially effective doses (min(X))
                                    LTY_A1[i,j,k,d1,d2] = LTY_int(X) + Yield[i,j,k,d1,d2] # LTY from point we arrive at/ Yield contribution from this year. Bellman step
                                    print(LTY_int(X),Yield[i,j,k,d1,d2],f_no)
                                if f_no == 1 and Yield[i,j,k,d1,d2]>params.Yield_threshold: # first year just get max yield.
                                    LTY_A1[i,j,k,d1,d2] = Yield[i,j,k,d1,d2]
                        ##                        
                        # update LTY_A2/dose_array_LTY
                        if np.amax(LTY_A1[i,j,k,:,:])>0:
                            LTY_A2[i,j,k] = np.amax(LTY_A1[i,j,k,:,:])
                            index = np.argwhere((LTY_A1[i,j,k,:,:] == np.amax(LTY_A1[i,j,k,:,:])))
                            dose_array_LTY[i,j,k,0] = 0.5*index[0][0]/(n_d-1)
                            dose_array_LTY[i,j,k,1] = 0.5*index[0][1]/(n_d-1)
                        ##
    ##
    LTY_int_new   = RegularGridInterpolator((x,y,z),LTY_A2)
    dictionary = {'F': F_array, 'phi_vec': phi_vec, 'LTY_int_new': LTY_int_new,'LTY_A1': LTY_A1,'LTY_A2': LTY_A2, 'D_A_LTY': dose_array_LTY}
    return dictionary
#----------------------------------------------------------------------------------------------
if __name__ == "__main__":
    asex_dict = {'n_phi': 3,
    'n_d':       4,
    'n_recurs': 15,
    'cube_arg':  0,
    'save_R0':   0,
    'plots':     1
    }
    asex_dict['load_R0'] = 1 - asex_dict['save_R0']
    asex_dict['save_F']  = asex_dict['save_R0']
    asex_dict['load_F']  = asex_dict['load_R0']

    #----------------------------------------------------------------------------------------------
    directory_string = 'Nick/HR_HR/Asexual_output/Saved_pickles/Old/'
    param_string   = ',n_phi=' + str(asex_dict['n_phi'])
    param_string_F = ',n_phi=' + str(asex_dict['n_phi']) + ',n_d=' + str(asex_dict['n_d']) + ',n_R=' + str(asex_dict['n_recurs']) + ',cube_arg=' + str(asex_dict['cube_arg'])
    R0_string = directory_string + 'R0_logged' + param_string + '.pickle'
    F_string  = directory_string + 'F_logged' + param_string_F + '.pickle'
    #----------------------------------------------------------------------------------------------
    if asex_dict['save_R0'] == 1:
        output  = FD_space(n_p=asex_dict['n_phi'])
        F       = output['F_array']
        LTY_int = output['Int_Y2'] # full dose gives max LTY for final year
        with open(R0_string, 'wb') as handle:
            pickle.dump(output,handle,protocol=pickle.HIGHEST_PROTOCOL) # protocol?
    if asex_dict['load_R0'] == 1:
        with open(R0_string, 'rb') as handle:
            output = pickle.load(handle)
            F        = output['F_array']
            LTY_int  = output['Int_Y2'] # full dose gives max LTY for final year
            Y2_array = output['Y2_array']
            Interpolate_dictionary = {**output, **params_dict}


    #----------------------------------------------------------------------------------------------
    # Have F_0
    if asex_dict['save_F'] == 1:
        dose_array_LTY = -1*np.ones((asex_dict['n_phi'],asex_dict['n_phi'],asex_dict['n_phi'],2))
        LTY_A2 = np.zeros((asex_dict['n_phi'],asex_dict['n_phi'],asex_dict['n_phi']))
        for i in range(1,asex_dict['n_recurs']):
            F_dict = recursion(F,f_no=i,dose_array_LTY = dose_array_LTY,LTY_A2=LTY_A2,LTY_int=LTY_int,n_d=asex_dict['n_d'],cube_arg=asex_dict['cube_arg'])
            dose_array_LTY = F_dict['D_A_LTY']
            LTY_A2         = F_dict['LTY_A2']
            LTY_int        = F_dict['LTY_int_new']
            print(LTY_A2)
            # print(LTY_int([-2,-2,-2]))
        ##
        F_final = F_dict['F']
        phi_vec = F_dict['phi_vec']
        LTY_A2  = F_dict['LTY_A2']
        D_A_LTY = F_dict['D_A_LTY']
        ##
        with open(F_string, 'wb') as handle:
            pickle.dump(F_dict,handle,protocol=pickle.HIGHEST_PROTOCOL) # protocol?
    if asex_dict['load_F'] == 1:
        with open(F_string, 'rb') as handle:
            F_load = pickle.load(handle)
            F_final = F_load['F']
            phi_vec = F_load['phi_vec']
            LTY_A2  = F_load['LTY_A2']
            D_A_LTY = F_load['D_A_LTY']
            F_dictionary = {**F_load, **params_dict}

#----------------------------------------------------------------------------------------------
    if asex_dict['plots'] == 1:
        # issues here
        cube_scatter_plot(F_final)
        ##-------------------------------------------
        mosaic_plot_2d(Array_to_plot=LTY_A2,inline_vec_to_plot =inline_vec_LTY,fig_type='con',title='LTY_A2')
        mosaic_plot_2d(Array_to_plot=F_final,inline_vec_to_plot=inline_vec,fig_type='con',title='F_final')
        ##---------------------------------------------------------------------------           
        mosaic_plot_2d(D_A_LTY,index=0,title='Dose_F0')
        mosaic_plot_2d(D_A_LTY,index=1,title='Dose_F1')
        # #----------------------------------------------------------------------------------------------
        # for i in range(asex_dict['n_phi']):
        #     Overlay_plotter(Con_mats=(F_final[:,i,:]),Con_levels=(inline_vec),Con_inline=('inline'),title=r'$log_{10} (p_{sr}) = $' + '%s' % phi_vec[i], x_lab=r'$log_{10} (p_{rs})$',y_lab=r'$log_{10} (p_{rr})$',xtick=phi_vec,ytick=phi_vec)
        #     # print(F_final[:,i,:])
        #----------------------------------------------------------------------------------------------            
        Z = generate_array(F_final)
        surface_plot_3d(Z)
        ###
        surface_plot_3d(Z,multiplot=0)
        #----------------------------------------------------------------------------------------------
        plt.show()