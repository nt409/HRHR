import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from math import ceil, log10
import copy
# import pickle
# import json
# import os
##
# import sys
# sys.path.insert(0, '/utils/')
##
from utils.Optimal_simulator_functions_and_animator import FD_space
from utils.functions_HRHR import object_dump, object_open
from utils.parameters_HRHR import params, params_dict
from utils.plotter_HRHR import mosaic_plot_2d, cube_scatter_plot, surface_plot_3d, generate_array, dose_choice_plotter, phi_space_navigator
from run2_optimal_asexual_cluster import param_string_rec
#----------------------------------------------------------------------------------------------
a_s_dict = object_open(params.JSON_path+'global.json','json')
##
a_s_dict['phi_vec']   = np.linspace(a_s_dict['log_phi_min'],0,a_s_dict['n_phi'])
a_s_dict['phi_vec_rr']   = np.linspace(2*a_s_dict['log_phi_min'],0,2*a_s_dict['n_phi'])
#----------------------------------------------------------------------------------------------
param_string_R0 = ',n_phi=' + str(a_s_dict['n_phi']) + ',log_phi_min=' + str(a_s_dict['log_phi_min'])
R0_string   = params.pickle_path + 'R0_logged' + param_string_R0 + '.pickle'
comb_string = params.pickle_path + 'Combined' + param_string_rec + '.pickle'
Final_output_string = params.pickle_path + 'Final_output' + param_string_rec + '.pickle'
#----------------------------------------------------------------------------------------------
if a_s_dict['load_comb'] != 1:
    Yield_out, prr2_out, prs2_out, psr2_out = [np.zeros((2*a_s_dict['n_phi'],a_s_dict['n_phi'],a_s_dict['n_phi'],a_s_dict['n_d'],a_s_dict['n_d'])) for j in range(4)]
    ##
    for k in range(2*a_s_dict['n_phi']):
        a_s_dict['phi_rr_val'] = a_s_dict['phi_vec_rr'][k]
        rec_string  = params.pickle_path + 'rec_logged' + param_string_rec + ',phi_rr_val=' + str(round(a_s_dict['phi_rr_val'],2)) + '.pickle'
        rec_load = object_open(rec_string)
        Yield_out[k,:,:,:,:] = rec_load['Yield']
        prr2_out[k,:,:,:,:]  = rec_load['prr2']
        prs2_out[k,:,:,:,:]  = rec_load['prs2']
        psr2_out[k,:,:,:,:]  = rec_load['psr2']
    ##
    comb_output = {'Yield_out': Yield_out,'prr2_out': prr2_out,'prs2_out': prs2_out,'psr2_out': psr2_out}
    
    object_dump(comb_string,comb_output)
else:
    comb_output = object_open(comb_string)

## print parameters
# a_s_dict['phi_rr_val'] = a_s_dict['phi_vec_rr'][3]
# rec_string  = params.pickle_path + 'rec_logged' + param_string_rec + ',phi_rr_val=' + str(round(a_s_dict['phi_rr_val'],2)) + '.pickle'
# rec_load2 = object_open(rec_string)
# print(rec_load2)

#----------------------------------------------------------------------------------------------
if a_s_dict['load_R0'] != 1:
    R0_output  = FD_space(comb_output['Yield_out'],bottom_phi=a_s_dict['log_phi_min'])
    F       = R0_output['F_array']
    LTY_int = R0_output['Int_Y2'] # full dose gives max LTY for final year
    
    object_dump(R0_string,R0_output)
else:
    R0_output = object_open(R0_string)
    F        = R0_output['F_array']
    LTY_int  = R0_output['Int_Y2'] # full dose gives max LTY for final year
    Y2_array = R0_output['Y2_array']
    Interpolate_dictionary = {**R0_output, **a_s_dict, **params_dict}


#----------------------------------------------------------------------------------------------
def recursion_algorithm(prr2,prs2,psr2,Yield_out,F_array,f_no,LTY_opt,LTY_int,dose_array_LTY,problem_points,log_phi_min = a_s_dict['log_phi_min']):
    print(f_no)
    n_p = prr2.shape[1]
    n_d = prr2.shape[-1]
    phi_vec = np.linspace(log_phi_min,0,n_p)
    x     = np.linspace(2*log_phi_min,0,2*n_p)
    y,z   = [phi_vec for i in range(2)]
    Int_F   = RegularGridInterpolator((x,y,z), F_array, bounds_error=False, fill_value=f_no)
    ##
    LTY_all      = np.zeros((2*n_p,n_p,n_p,n_d,n_d))
    LTY_opt      = copy.deepcopy(LTY_opt) # interpolation function was changing when the array did
    for i in range(2*n_p):
        for j in range(n_p):
            for k in range(n_p):
                for d1 in range(n_d):
                    for d2 in range(n_d):
                        X = [log10(prr2[i,j,k,d1,d2]),log10(prs2[i,j,k,d1,d2]),log10(psr2[i,j,k,d1,d2])]
                        if Yield_out[i,j,k,d1,d2]>params.Yield_threshold: # so that haven't already been 'upgraded'
                            if X[0]<2*log_phi_min or X[0]>0:
                                problem_points[i,j,k,d1,d2,0] = 0 # if the point gets sent out of the cube on side iii
                                print(X)
                                print(i,j,k,f_no)
                            for iii in [1,2]: # if out of range, project on to nearest point on cube
                                if X[iii]<log_phi_min or X[iii]>0:
                                    problem_points[i,j,k,d1,d2,iii] = 0 # if the point gets sent out of the cube on side iii
                                    print(X)
                                    print(i,j,k,f_no)
                                    # project back onto cube?
                                    if X[iii]<log_phi_min:
                                        X[iii]=log_phi_min
                                    if X[iii]>0:
                                        X[iii]=0

                        if Yield_out[i,j,k,d1,d2]>params.Yield_threshold and F_array[i,j,k]==f_no: # so that haven't already been 'upgraded'                        
                            proceed = True
                            if problem_points[i,j,k,d1,d2,iii] == 0: # if the point gets sent out of the cube on side iii
                                proceed = False
                            if proceed:
                                # update F
                                if Int_F(X) > f_no - 0.01: # if we get sent to somewhere with at least as good a f_no # if min(X)>log_phi_min and max(X)<0 and
                                    F_array[i,j,k] = f_no + 1

                                # update LTY_all
                                # first year just get max yield.
                                if f_no == 1:
                                    LTY_all[i,j,k,d1,d2] = Yield_out[i,j,k,d1,d2]
                                # subsequent years
                                if F_array[i,j,k]>=f_no and f_no>1: # if we can survive at least f_no years # if min(X)>log_phi_min and max(X)<0
                                    LTY_all[i,j,k,d1,d2] = LTY_int(X) + Yield_out[i,j,k,d1,d2] # LTY from point we arrive at/ Yield contribution from this year. Bellman step
                        ##
                # update LTY_opt/dose_array_opt_LTY
                if np.amax(LTY_all[i,j,k,:,:])>0:
                    LTY_opt[i,j,k] = np.amax(LTY_all[i,j,k,:,:])
                    index = np.argwhere((LTY_all[i,j,k,:,:] == np.amax(LTY_all[i,j,k,:,:])))
                    dose_array_LTY[i,j,k,0] = 0.5*index[0][0]/(n_d-1)
                    dose_array_LTY[i,j,k,1] = 0.5*index[0][1]/(n_d-1)
    LTY_int_new   = RegularGridInterpolator((x,y,z),LTY_opt)
    dictionary = {'F': F_array, 'phi_vec': phi_vec, 'LTY_int_new': LTY_int_new,'LTY_all': LTY_all,'LTY_opt': LTY_opt, 'D_A_LTY': dose_array_LTY, 'P_points': problem_points}
    return dictionary

#----------------------------------------------------------------------------------------------
n_p = a_s_dict['n_phi']
n_d = a_s_dict['n_d']
# Have F_0
if a_s_dict['load_final']==0:
    dose_array_LTY = -1*np.ones((2*n_p,n_p,n_p,2))
    problem_points = np.ones((2*n_p,n_p,n_p,n_d,n_d,3))
    LTY_opt = np.zeros((2*n_p,n_p,n_p))
    for f_no in range(1,a_s_dict['n_recurs']):
        F_dict = recursion_algorithm(comb_output['prr2_out'],comb_output['prs2_out'],comb_output['psr2_out'],comb_output['Yield_out'],F,f_no,LTY_opt,LTY_int,dose_array_LTY,problem_points)
        dose_array_LTY = F_dict['D_A_LTY']
        LTY_opt        = F_dict['LTY_opt']
        LTY_int        = F_dict['LTY_int_new']
        F              = F_dict['F']
        problem_points = F_dict['P_points']
        plt.show()
    object_dump(Final_output_string,F_dict)
if a_s_dict['load_final']==1:
    F_dict = object_open(Final_output_string,'pickle')


x     = a_s_dict['phi_vec_rr']
y,z   = [a_s_dict['phi_vec'] for i in range(2)]
Dose_int_1 = RegularGridInterpolator((x,y,z),F_dict['D_A_LTY'][:,:,:,0])
Dose_int_2 = RegularGridInterpolator((x,y,z),F_dict['D_A_LTY'][:,:,:,1])
F_int = RegularGridInterpolator((x,y,z),F_dict['F'])


if a_s_dict['plots'] == 1:
    cube_plot = True
    dose_plot = False
    mosaics = False
    plot_3d = False
    
    inline_vec = np.concatenate((range(a_s_dict['n_recurs']+1),[a_s_dict['n_recurs']+0.999]))
    inline_vec_LTY = np.linspace(0,100*a_s_dict['n_recurs'],a_s_dict['n_recurs']+1)
    inline_vec_DA = np.concatenate(([-0.1],np.linspace(0,0.5,101)))
    ##
    L,R,B,T = (0.12,0.9,0.12,0.9)
    marker_index = ['o','<','>','v','^','s','d','x','+','.',',','o','<','>','v','^','s','d','x','+','.',',']
    alpha = 0.6
    # print(F_dict['F'])
    if cube_plot:
        cube_scatter_plot(F_dict['F'],phi_vec_to_plot=a_s_dict['phi_vec'],phi_rr_vec=a_s_dict['phi_vec_rr'])
    ##-------------------------------------------
    if dose_plot:
        # problem?
        Z0 = [-4,-4,-4]
        n_years = ceil(F_int(Z0))
        if n_years<1:
            n_years = a_s_dict['n_recurs']


        ZZ = np.zeros((3,n_years))
        print(ZZ,Z0)
        ZZ[:,0] = Z0
        doses = -1*np.ones((2,n_years-1))
        Yield = np.zeros(n_years-1)

        for i in range(n_years-1):
            ZZ[:,i+1], doses[:,i], Yield[i] = phi_space_navigator(ZZ[:,i],Dose_int_1,Dose_int_2)
        dose_choice_plotter(ZZ,doses,Yield,a_s_dict['phi_vec'],a_s_dict['phi_vec_rr'])
    ##-------------------------------------------
    probs = np.ones((2*n_p,n_p,n_p,3))
    for i in range(2*n_p):
        for j in range(n_p):
            for k in range(n_p):
                probs[i,j,k,0] = np.amin(F_dict['P_points'][i,j,k,:,:,0])
                probs[i,j,k,1] = np.amin(F_dict['P_points'][i,j,k,:,:,1])
                probs[i,j,k,2] = np.amin(F_dict['P_points'][i,j,k,:,:,2])
    if mosaics:
        mosaic_plot_2d(Array_to_plot=F_dict['LTY_opt'],phi_vec_to_plot=a_s_dict['phi_vec'],phi_rr_vec=a_s_dict['phi_vec_rr'],inline_vec_to_plot =inline_vec_LTY,fig_type='con',title='LTY_opt')
        mosaic_plot_2d(Array_to_plot=probs[:,:,:,0],phi_vec_to_plot=a_s_dict['phi_vec'],phi_rr_vec=a_s_dict['phi_vec_rr'],inline_vec_to_plot =[0,0.5,1],fig_type='col',title='P_points_rr',colorbar_on=False)
        mosaic_plot_2d(Array_to_plot=probs[:,:,:,1],phi_vec_to_plot=a_s_dict['phi_vec'],phi_rr_vec=a_s_dict['phi_vec_rr'],inline_vec_to_plot =[0,0.5,1],fig_type='col',title='P_points_rs',colorbar_on=False)
        mosaic_plot_2d(Array_to_plot=probs[:,:,:,2],phi_vec_to_plot=a_s_dict['phi_vec'],phi_rr_vec=a_s_dict['phi_vec_rr'],inline_vec_to_plot =[0,0.5,1],fig_type='col',title='P_points_sr',colorbar_on=False)
        mosaic_plot_2d(Array_to_plot=F_dict['F']   ,phi_vec_to_plot=a_s_dict['phi_vec'],phi_rr_vec=a_s_dict['phi_vec_rr'],inline_vec_to_plot=inline_vec,fig_type='con',title='F_final')
        ##---------------------------------------------------------------------------           
        ind_v =  [0,50,60,64,68,72,77,80,86]

        mosaic_plot_2d(Array_to_plot=F_dict['D_A_LTY'],phi_vec_to_plot=a_s_dict['phi_vec'],phi_rr_vec=a_s_dict['phi_vec_rr'],inline_vec_to_plot=inline_vec_DA,fig_type='col',index=0,title='Optimal dose of fungicide 1',index_vec=ind_v)
        inline_vec_LTY = [np.amin(F_dict['LTY_opt']),np.amax(F_dict['LTY_opt'])]
        inline_vec     = [np.amin(F_dict['F']),np.amax(F_dict['F'])]
        mosaic_plot_2d(Array_to_plot=F_dict['LTY_opt'],phi_vec_to_plot=a_s_dict['phi_vec'],phi_rr_vec=a_s_dict['phi_vec_rr'],inline_vec_to_plot =inline_vec_LTY,fig_type='col',title='Optimal Lifetime Yield',index_vec=ind_v)
        mosaic_plot_2d(Array_to_plot=F_dict['F']   ,phi_vec_to_plot=a_s_dict['phi_vec'],phi_rr_vec=a_s_dict['phi_vec_rr'],inline_vec_to_plot=inline_vec,fig_type='col',title='Optimal Effective Life',index_vec=ind_v)



        # mosaic_plot_2d(Array_to_plot=F_dict['D_A_LTY'],phi_vec_to_plot=a_s_dict['phi_vec'],phi_rr_vec=a_s_dict['phi_vec_rr'],inline_vec_to_plot=inline_vec_DA,fig_type='col',index=0,title='Optimal dose of fungicide 1',index_vec=ind_v2)
        # mosaic_plot_2d(Array_to_plot=F_dict['D_A_LTY'],phi_vec_to_plot=a_s_dict['phi_vec'],phi_rr_vec=a_s_dict['phi_vec_rr'],inline_vec_to_plot=inline_vec_DA,fig_type='col',index=1,title='Optimal dose of fungicide 2')
    # #----------------------------------------------------------------------------------------------
    # for i in range(a_s_dict['n_phi']):
    #     Overlay_plotter(Con_mats=(F_final[:,i,:]),Con_levels=(inline_vec),Con_inline=('inline'),title=r'$log_{10} (p_{sr}) = $' + '%s' % a_s_dict['phi_vec'][i], x_lab=r'$log_{10} (p_{rs})$',y_lab=r'$log_{10} (p_{rr})$',xtick=a_s_dict['phi_vec'],ytick=a_s_dict['phi_vec'])
    #     # print(F_final[:,i,:])
    #----------------------------------------------------------------------------------------------            
    if plot_3d:
        Z = generate_array(F_dict['F'],phi_vec_here=a_s_dict['phi_vec'],phi_rr_vec=a_s_dict['phi_vec_rr'])
        surface_plot_3d(Z,phi_vec_to_plot=a_s_dict['phi_vec'],phi_rr_vec=a_s_dict['phi_vec_rr'])
        # ###
        surface_plot_3d(Z,phi_vec_to_plot=a_s_dict['phi_vec'],phi_rr_vec=a_s_dict['phi_vec_rr'],multiplot=False)

    plt.show()