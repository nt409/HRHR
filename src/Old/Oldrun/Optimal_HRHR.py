import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pickle
from math import ceil, log, floor, log10, exp
##
# import sys
# sys.path.insert(0, '/utils/')
##
import os
##
from utils.params import params, params_dict
from utils.functions import interpolate, master_loop_one_tactic, master_loop_grid_of_tactics
# from utils.plotter_HRHR import Overlay_plotter, Season_plotter, Equal_dose, plot_disease_dynamics, bar_plotter
from utils.Optimal_simulator_functions_and_animator import optimal_simulator, optimal_animator, optimal_mosaic, phi_plane_jumps, transformed_yield_contours, Y_t_space_contour_finder, Int_Y_FD, C_contours
#----------------------------------------------------------------------------------------------
##
p    = 0.01 # keep in (0,0.5)
phi  = 2*10**(-2)
##
res1 = p*phi # r1 = 10**(-3) # r2 = 10**(-1)
res2 = (1-p)*phi
num_phi = 8
n_p_i   = 2*num_phi
##
n_d   = 5 # 8      # n_doses
n_d_i = 2*n_d # 2*n_d   # n_doses_interpolated # maybe limit to 2*n_doses?
n_seasons = 9
##
interp                = True
equal_dose_on_contour = True
int_type              = 'linear' # should only use cubic with enough data, ideally around 8 plus
g_or_y                = 'G' # 'D', 'Y', 'G'
#----------------------------------------------------------------------------------------------
opt_HRHR_params = {'opt_p': p,'opt_phi': phi,'opt_res1': res1,'opt_res2': res2,'opt_n_d': n_d,'opt_n_d_i': n_d_i,'opt_n_seasons': n_seasons,'opt_interp': interp,'opt_equal_dose_on_contour': equal_dose_on_contour,'opt_int_type': int_type,'opt_g_or_y': g_or_y}
#----------------------------------------------------------------------------------------------
directory_string = 'C:/Users/user/Documents/Python/Nick/HR_HR/Asexual_output/Saved_pickles/Old/'
##
param_string_I_Y_FD = ',nphi='+str(num_phi)+',n_int='+str(n_p_i)+'.pickle'
param_string_cdg = ',nd='+str(n_d)+',p='+str(p)+',phi='+str(phi)+'.pickle'
param_string_opt_s = ',nd='+str(n_d)+',ndi='+str(n_d_i)+',p='+str(p)+',phi='+str(phi)+'.pickle'
##
I_Y_FD_string = directory_string+'I_Y_FD'+param_string_I_Y_FD
cdg_string    = directory_string+'con_dose_grid'+param_string_cdg
opt_s_string  = directory_string+'opt_s'+param_string_opt_s
##
#----------------------------------------------------------------------------------------------
## Modify these
##
phi_cont_plt     = False   # poster fig 3 # plot with phi_1 phi_2 values for diff seasons and 95 contour for full dose.
Y_diff_plot      = False   # poster fig 4
c_cont_plt       = False   # horizontal thing as phi_1 phi_2 varies # needs run_optimal_sim
phi_cont_plt_2   = False   # plot showing where we move from a particular value of phi_1, phi_2 (either all options or smallest number of contour jumps)
phi_plane_jump_2 = False   # overlay_plotter problem   # can run on its own # plot showing where we move from particular value, by colour
pp_plots         = False   # overlay_plotter problem
Y_t_space        = False   # No_R_pt problem
##
present_plots    = False  # as in FY presentation
season_plt       = False # Overlay_plotter problem  # plots RF_ratio, Y1, Y1_int by season   
contour_plt      = False # Overlay_plotter problem  # plots individual seasons by G  
overlay_plt      = False # Overlay_plotter problem  # plots individual seasons by G  
surface_plt      = False                            # plots individual seasons by G  
mosaic_plt       = True # poster fig 2 # plots in dose space by year, shows selected tactic
##
animate          = False
##
#----------------------------------------------------------------------------------------------
## # shouldn't need to modify these
plots                  = False         # turn all plots here off # needs run_optimal_sim
phi_plane_jump_1       = False    # can run on its own - needed if using diff
run_interp_Y_FD_by_phi = False # needed for Ytspace, opt_sim
run_const_grid         = False # needed for phi_cont_plt_2
run_optimal_sim        = False # run_optimal_sim needed for all below. Needs phi_plane_jump_1 for Diff, run_interp_Y_FD_by_phi
##
if present_plots or season_plt or contour_plt or overlay_plt or surface_plt or mosaic_plt:
    plots = True
if phi_cont_plt_2:
    run_const_grid = True
if g_or_y == 'D' or Y_diff_plot or mosaic_plt:
    phi_plane_jump_1 = True
if plots or phi_cont_plt or c_cont_plt or animate: # present_plots or season_plt or contour_plt or overlay_plt or surface_plt or mosaic_plt
    run_optimal_sim = True
if Y_t_space or run_optimal_sim:
    run_interp_Y_FD_by_phi = True
##
# #----------------------------------------------------------------------------------------------
##
if run_interp_Y_FD_by_phi:
    if os.path.exists(I_Y_FD_string):
        with open(I_Y_FD_string, 'rb') as handle:
            loaded_int = pickle.load(handle)
            int_dictionary = {**loaded_int, **params_dict, **opt_HRHR_params}
            Interpolate_Y_at_FD = loaded_int['Int_Y_FD']
    else:
        Interpolate_Y_at_FD, Int_Y_at_FD = Int_Y_FD(n_phi=num_phi,output_size=n_p_i)
        interp_Y_FD_dict = {'Int_Y_FD': Interpolate_Y_at_FD}
        with open(I_Y_FD_string, 'wb') as handle:
            pickle.dump(interp_Y_FD_dict,handle,protocol=pickle.HIGHEST_PROTOCOL) # protocol?
# #----------------------------------------------------------------------------------------------
if run_const_grid:
    if os.path.exists(cdg_string):
        with open(cdg_string, 'rb') as handle:
            const_dose_grid = pickle.load(handle)
            const_dose_grid = {**const_dose_grid, **params_dict, **opt_HRHR_params}
    else:
        const_dose_grid = master_loop_grid_of_tactics(n_d,n_seasons,res1,res2)
        with open(cdg_string, 'wb') as handle:
            pickle.dump(const_dose_grid,handle,protocol=pickle.HIGHEST_PROTOCOL) # protocol?
    
    print(const_dose_grid['FY'])
# #----------------------------------------------------------------------------------------------
if phi_cont_plt_2: # needs const_dose_grid
    n = n_d
    dose = 0.5
    C = 0.1 # 0.13
    ##
    YieldV_phi_cont_plt2 = np.zeros((n,n))
    xy_array = np.zeros((n,n))
    phi_1,phi_2     = [np.linspace(0,1,n) for i in range(2)]
    Rv1, Rv2        = [-1*np.ones(n_d**2 + 1) for i in range(2)]
    Rv1_use,Rv2_use = [-1*np.ones(2) for i in range(2)]
    ##
    for i in range(len(phi_1)):
        for j in range(len(phi_2)):
            YieldV_phi_cont_plt2[i,j] = master_loop_one_tactic([dose],[dose],[dose],[dose],phi_1[i],phi_2[j])['Yield_vec']
    output_size = n_d
    Interpolate_Y, Int_Y = interpolate(YieldV_phi_cont_plt2,len(YieldV_phi_cont_plt2),output_size)
    
    top = 0
    Rv1[0] = res1
    Rv2[0] = res2
    Rv1_use[0] = res1
    Rv2_use[0] = res2
    k = 1
    for i in range(n_d):
        for j in range(n_d):
            if const_dose_grid['Yield'][i,j,1]>params.Yield_threshold:
                Rv1[k] = const_dose_grid['Res_array_1'][i,j,1]
                Rv2[k] = const_dose_grid['Res_array_2'][i,j,1]
                Y2 = Interpolate_Y(Rv1[k],Rv2[k])
                k = k + 1
                if Y2 > top:
                    top = Y2

    for i in range(n_d):
        for j in range(n_d):
            if top == Interpolate_Y(const_dose_grid['Res_array_1'][i,j,1],const_dose_grid['Res_array_2'][i,j,1]):
                Rv1_use[1] = const_dose_grid['Res_array_1'][i,j,1]
                Rv2_use[1] = const_dose_grid['Res_array_2'][i,j,1]         
    ##
    # cont = [params.Yield_threshold,99,99.05,99.06,99.0668,99.07]
    cont = [params.Yield_threshold]

    Con_attrs = {'Con_mats': (YieldV_phi_cont_plt2),'Con_levels':(cont),'Con_inline':('inline')}
    Scat_att  = {'x_scat':Rv1,'y_scat':Rv2}
    Overlay_plotter(Con_attr=Con_attrs,figure_attributes={'x_lab':r'$\phi_1$','y_lab':r'$\phi_2$'},Scat_attr= Scat_att)#,title=r'Yield, $\phi_1$, $\phi_2$'

    Con_attrs2 = {'Con_mats': (YieldV_phi_cont_plt2),'Con_levels':(cont),'Con_inline':('inline')}
    Scat_att2  = {'x_scat':Rv1_use,'y_scat':Rv2_use}
    Overlay_plotter(Con_attr=Con_attrs2,figure_attributes={'x_lab':r'$\phi_1$','y_lab':r'$\phi_2$','title':r'Yield, $\phi_1$, $\phi_2$'},Scat_attr= Scat_att2)#,title=r'Yield, $\phi_1$, $\phi_2$'

# #----------------------------------------------------------------------------------------------
if phi_plane_jump_1:
    grid    = False
    n_doses = 7 # 15 #4 #6 # number of doses to iterate over
    ##
    n_phi   = 3 #4 #10 # if grid = 1, number of phi values to iterate over 
    n_steps = 2 # if grid = 0
    ##
    dose_descriptors = False # if 1, draw diamonds with colours corresponding to full/zero dose
    phi_vector = [1] # range(n_steps) # range(n_steps) # [0] # [n_steps-1]
    c_cont  = np.linspace(98,99.07,3) # [99]
    ##
    param_string_pp1 = ',nd='+str(n_doses)+',n_phi='+str(n_phi)+',n_steps='+str(n_steps)+'.pickle'
    pp1_string = directory_string+'pp1'+param_string_pp1
    ##
    if os.path.exists(pp1_string):
        with open(pp1_string, 'rb') as handle:
            pp_1 = pickle.load(handle)
            pp_1 = {**pp_1, **params_dict, **opt_HRHR_params}
    else:
        pp_1 = phi_plane_jumps(phi_grid=1,n_phi=n_phi,n_d=n_doses,n_steps=n_steps)
        with open(pp1_string, 'wb') as handle:
            pickle.dump(pp_1,handle,protocol=pickle.HIGHEST_PROTOCOL)
# #----------------------------------------------------------------------------------------------
    if Y_diff_plot:
        for i in c_cont:
            transformed_yield_contours(grid=grid,Y_cont=i,dose_descriptors=dose_descriptors,n_steps=n_steps,n_phi=n_phi,n_doses=n_doses,phi_vector=phi_vector,x_lim=(params.Yield_threshold-1,99.5),title=False)
        ##
        cont_yield = np.linspace(97.5,99,7)
        cont_yield2 = np.concatenate(([params.Yield_threshold],np.linspace(97.5,99,4),[100]))
        cont3 = np.concatenate((np.linspace(0,0.2,5),[0.5,1,2,3]))


        # Overlay_plotter(Con_mats=(pp_1['Diff'],pp_1['Yield_vec'],pp_1['Yield_vec']),Contourf=(None,'k',None),Con_inline=('inline',None,'Yield'),Con_levels=(cont3,[0,params.Yield_threshold],cont_yield),x_lab=r'$\phi_1$',y_lab=r'$\phi_2$')
        # Overlay_plotter(Con_mats=(pp_1['Yield_vec'],pp_1['Yield_vec'],pp_1['Diff']),Contourf=(1,'k',None),Con_inline=('Yield',None,'inline'),Con_levels=(cont_yield2,[0,params.Yield_threshold],cont3),x_lab=r'$\phi_1$',y_lab=r'$\phi_2$')
        # Overlay_plotter(Con_mats=(pp_1['Diff'],pp_1['Yield_vec']),Contourf=(None,'k'),Con_inline=('inline',None),Con_levels=(cont3,[0,params.Yield_threshold]),x_lab=r'$\phi_1$',y_lab=r'$\phi_2$')
        # Overlay_plotter(Con_mats=(pp_1['Diff'],pp_1['Yield_vec']),Contourf=(1,'k'),Con_inline=(None,None),Con_levels=(cont3,[0,params.Yield_threshold]),x_lab=r'$\phi_1$',y_lab=r'$\phi_2$')
        # Overlay_plotter(Con_mats=(pp_1['Diff'],pp_1['Yield_vec'],pp_1['Yield_vec']),Contourf=(1,'k',None),Con_inline=('Difference',None,'inline'),Con_levels=(cont3,[0,params.Yield_threshold],cont_yield),x_lab=r'$\phi_1$',y_lab=r'$\phi_2$')
# #----------------------------------------------------------------------------------------------
if phi_plane_jump_2:
    n_doses = 8 # 10 # number of doses to iterate over
    n_phi   = 10 # 14 # 8 # number of phi values to iterate over
    n_steps = 3 # 6 # 4 # number of points along contour to pick
    C_cont  = 99.071 # 99.05 # FD Yield contour to consider
    ##
    param_string_pp2 = ',nd='+str(n_doses)+',n_phi='+str(n_phi)+',n_steps='+str(n_steps)+',C_cont='+str(C_cont)+'.pickle'    
    pp2_string = directory_string+'pp2'+param_string_pp2
    ##
    if os.path.exists(pp2_string):
        with open(pp2_string, 'rb') as handle:
            pp_2 = pickle.load(handle)
            pp_2 = {**pp_2, **params_dict, **opt_HRHR_params}
    else:
        pp_2 = phi_plane_jumps(C_cont=C_cont,n_phi=n_phi,n_d=n_doses,n_steps=n_steps)
        with open(pp2_string, 'wb') as handle:
            pickle.dump(pp_2,handle,protocol=pickle.HIGHEST_PROTOCOL)
    ##
    if pp_plots:
        for k in range(1,n_steps):
            Overlay_plotter(Con_mats=(pp_2['SR1'][:,:,k],pp_2['Yield'][:,:,0,k]),Con_inline=('inline',None),Contourf=(None,'k'),Con_levels=(None,[0,params.Yield_threshold]),title=r'SR1; $\phi_1$ %s, $\phi_2$ %s' % (pp_2['phi_1_values'][k],pp_2['phi_2_values'][k]))
            Overlay_plotter(Con_mats=(pp_2['SR2'][:,:,k],pp_2['Yield'][:,:,0,k]),Con_inline=('inline',None),Contourf=(None,'k'),Con_levels=(None,[0,params.Yield_threshold]),title=r'SR2; $\phi_1$ %s, $\phi_2$ %s' % (pp_2['phi_1_values'][k],pp_2['phi_2_values'][k]))
        for k in range(2,n_steps):
            # Overlay_plotter(Con_mats=(pp_2['dist1'][:,:,k],pp_2['Yield'][:,:,0,k]),Con_inline=('inline',None),Contourf=(None,'k'),Con_levels=(np.linspace(0,np.amax(pp_2['dist1'][:,:,k]),5),[0,params.Yield_threshold]),title=r'Dist1; $\phi_1$ %s, $\phi_2$ %s' % (pp_2['phi_1_values'][k],pp_2['phi_2_values'][k]))
            # Overlay_plotter(Con_mats=(pp_2['dist1'][:,:,k],pp_2['Yield'][:,:,0,k]),Con_inline=('inline',None),Contourf=(None,'k'),Con_levels=(np.linspace(0,np.amax(pp_2['dist1'][:,:,k]),5),[0,params.Yield_threshold]),title=r'Dist1; $\phi_1$ %s, $\phi_2$ %s' % (pp_2['phi_1_values'][k],pp_2['phi_2_values'][k]))
            contz = np.linspace(0,max(np.amax(pp_2['dist1'][:,:,k]),np.amax(pp_2['dist2'][:,:,k])),5)
            Overlay_plotter(Con_mats=(pp_2['dist1'][:,:,k],pp_2['dist2'][:,:,k],pp_2['Yield'][:,:,0,k]),Con_inline=('inline',None,None),Contourf=(None,None,'k'),Con_levels=(contz,contz,[0,params.Yield_threshold]),title=r'Dist1/2; $\phi_1$ %s, $\phi_2$ %s' % (pp_2['phi_1_values'][k],pp_2['phi_2_values'][k]))
        for k in range(2,n_steps):
            Overlay_plotter(Con_mats=(pp_2['Y2'][:,:,k],pp_2['Yield'][:,:,0,k]),Con_inline=('inline',None),Contourf=(None,'k'),Con_levels=(np.linspace(params.Yield_threshold,max(ceil(np.amax(pp_2['Y2'][:,:,k])),96),11),[0,params.Yield_threshold]) ,title=r'Y2; $\phi_1$ %s, $\phi_2$ %s' % (pp_2['phi_1_values'][k],pp_2['phi_2_values'][k]))
    for kk in range(n_steps):
        print(pp_2['Yield'][:,:,0,kk],'... Yield')
        # Overlay_plotter(Con_mats=(pp_2['Y2'][:,:,kk]),Con_levels=(np.linspace(np.amax(pp_2['Y2'][:,:,kk])-1,np.amax(pp_2['Y2'][:,:,kk]),20)),Con_inline=('inline'),title='Y2 by dose, %s' % kk)
        Rv1_flat = pp_2['Rv1'][:,:,kk].flatten()
        Rv2_flat = pp_2['Rv2'][:,:,kk].flatten()
        cont =  pp_2['Y2'][:,:,kk].flatten()
        index = np.argwhere(cont==50)
        cont     = np.delete(cont,index)
        Rv1_flat = np.delete(Rv1_flat,index)
        Rv2_flat = np.delete(Rv2_flat,index)
        cont = np.sort(cont)
        C_min = cont[0]
        C_max = cont[-1]
        ##
        labels = [['Other','ZD_1','ZD_2','FD_1','FD_2','CD','FD_B'],['c','y','g','r','m','b','k']]
        cols = np.repeat([labels[1][0]],len(Rv1_flat))
        for k in range(len(Rv1_flat)):
            for j in range(n_doses):
                if pp_2['Rv1'][0,j,kk]==Rv1_flat[k] and j!=n_doses-1:
                    cols[k] = labels[1][1] # zd_1 colour
                if pp_2['Rv1'][j,0,kk]==Rv1_flat[k] and j!=n_doses-1:
                    cols[k] = labels[1][2] # zd_2 colour
                if pp_2['Rv1'][-1,j,kk]==Rv1_flat[k]:
                    cols[k] = labels[1][3] # fd_1 colour
                if pp_2['Rv1'][j,-1,kk]==Rv1_flat[k]:
                    cols[k] = labels[1][4] # fd_2 colour
                if j!=0 and j!=n_doses-1: # so not full dose of 2
                    v0= pp_2['Yield'][:,j,0,kk]
                    v = pp_2['Yield'][1:-1,j,0,kk] # so not full dose of 1
                    if min(v0)<params.Yield_threshold and max(v0)>params.Yield_threshold: # so the contour has a point for this j
                        v2 = np.delete(v,np.argwhere(v<=params.Yield_threshold)) # not 0 dose or full dose
                        for i in range(len(v0)):
                            if v0[i]==min(v2) and pp_2['Rv1'][i,j,kk]==Rv1_flat[k]: #cont[k] == pp_2['Y2'][i,j,kk]: # selects i on contour? Fixed dose of 2, dose 1 that gives us smallest yield
                                cols[k] = labels[1][5] # CD colour
            if pp_2['Rv1'][-1,-1,kk]==Rv1_flat[k]:
                cols[k] = labels[1][6] # fd_both colours
        ##
        cont = np.concatenate(([params.Yield_threshold],[C_min],[C_max]))
        cont2 = [params.Yield_threshold]
        phi_1_scat_k = [pp_2['phi_1_values'][kk]]
        phi_2_scat_k = [pp_2['phi_2_values'][kk]]
        Overlay_plotter(Con_mats=(pp_2['Yield_vec']),Con_inline=('Yield'),Con_levels=(cont),x_scat=Rv1_flat,y_scat=Rv2_flat,title=r'$p_1$, $p_2$ = %r, %r' % (round(pp_2['phi_1_values'][kk],4),round(pp_2['phi_2_values'][kk],4)),x_lab=r'$\phi_1$',y_lab=r'$\phi_2$',x_scat_k=phi_1_scat_k,y_scat_k=phi_2_scat_k,scat_colours=cols,scat_leg=labels) # 'Y_0,Y_M,Y_m,p_1,p_2 = %r, %r, %r, %r,%r' % (C_cont,round(C_max,2),round(C_min,2),round(pp_2['phi_1_values'][kk],2),round(pp_2['phi_2_values'][kk],2))
        # Overlay_plotter(Con_mats=(pp_2['Yield_vec']),Con_inline=('inline'),Con_levels=(cont2),x_scat=Rv1_flat,y_scat=Rv2_flat,title=r'Yield, $\phi_1$, $\phi_2$',x_lab=r'$\phi_1$',y_lab=r'$\phi_2$',x_scat_k=phi_1_scat_k,y_scat_k=phi_2_scat_k)
# #----------------------------------------------------------------------------------------------
if Y_t_space:
    n_season = 4 # 12
    dose_n   = 3
    t_vect = np.linspace(0,1,4)
    t_used, Y2_used, res_used, dose_used, bin_from = Y_t_space_contour_finder(Interpolate_Y_at_FD = Interpolate_Y_at_FD,t_vector=t_vect,n_seasons=n_season,n_d=dose_n,res1=res1,res2=res2)
    ##
    size = Y2_used.shape[0]*Y2_used.shape[1] + 1
    xp2, yp2, colour_vec, bin_f = np.zeros(size), np.zeros(size), np.zeros(size), np.zeros(size)
    xp2[0], yp2[0], colour_vec[0] = Interpolate_Y_at_FD(res1,res2), res1/(res1+res2), 0
    i = 1
    for j in range(n_season):
        for k in range(len(t_vect)-1):
            xp2[i]   = Y2_used[k,j]
            yp2[i]   = t_used[k,j]
            bin_f[i] = bin_from[k,j]
            colour_vec[i] = j
            i = i+1
    transformed_yield_contours(Interpolate_Y_FD=Interpolate_Y_at_FD,x_points=xp2,y_points=yp2,colour_vec=colour_vec,n_s=n_season,hlines=t_vect,bin_from=bin_f,log_ind=False)
# #----------------------------------------------------------------------------------------------
if run_optimal_sim != 1:
    plt.show()
# #----------------------------------------------------------------------------------------------
if run_optimal_sim:
    # defaults always_on_contour = 0, interp = 1, interp_res =1, interp_type = 'linear'... can change
    if os.path.exists(opt_s_string):
        with open(opt_s_string, 'rb') as handle:
            opt_s = pickle.load(handle)
            opt_s = {**opt_s, **params_dict, **opt_HRHR_params}
    else:
        opt_s = optimal_simulator(res1,res2,n_d,n_d_i,n_seasons,interp,Interpolate_Y_at_FD=Interpolate_Y_at_FD,interp_type= int_type,equal_dose_on_contour=equal_dose_on_contour,G_or_yield='D',Diff=pp_1['Diff'])
        with open(opt_s_string, 'wb') as handle:
            pickle.dump(opt_s,handle,protocol=pickle.HIGHEST_PROTOCOL)
    # optimal_sim_contour_output = optimal_simulator(res1,res2,n_d,n_d_i,n_seasons,interp,Interpolate_Y_at_FD=Interpolate_Y_at_FD,always_on_contour=1,G_or_yield=g_or_y,Diff=pp_1['Diff'])
    # this is to generate H,G,Y for colour plots, and the doses. Then we use these doses in a new simulation so no approximations
    # #----------------------------------------------------------------------------------------------
    dose_opt_not_int     = np.concatenate((opt_s['dose_array'][:,0:(opt_s['success_years'])],    0.5*np.ones((2,1))),axis=1) # plus one full dose to show final year impossible
    dose_opt_int = np.concatenate((opt_s['dose_array_int'][:,0:(opt_s['success_years'])],0.5*np.ones((2,1))),axis=1) # plus one full dose to show final year impossible
    ##
    d_const_full     = 0.5 *np.ones(opt_s['success_years']+1)
    d_const_med = 0.35*np.ones(opt_s['success_years']+1)
    d_const_low = 0.2 *np.ones(opt_s['success_years']+1)
    d_zero      =     np.zeros(opt_s['success_years']+1)
    # #----------------------------------------------------------------------------------------------
    not_int_1_tact    = master_loop_one_tactic(dose_opt_not_int[0,:],dose_opt_not_int[1,:],dose_opt_not_int[0,:],dose_opt_not_int[1,:],opt_s['res_array'][0,0],opt_s['res_array'][1,0])
    int_1_tact        = master_loop_one_tactic(dose_opt_int[0,:],dose_opt_int[1,:],dose_opt_int[0,:],dose_opt_int[1,:],opt_s['res_array'][0,0],opt_s['res_array'][1,0])
    ###
    const_1_tact_full = master_loop_one_tactic(d_const_full,d_const_full,d_const_full,d_const_full,opt_s['res_array'][0,0],opt_s['res_array'][1,0])
    const_1_tact_med  = master_loop_one_tactic(d_const_med,d_const_med,d_const_med,d_const_med,opt_s['res_array'][0,0],opt_s['res_array'][1,0])
    const_1_tact_low  = master_loop_one_tactic(d_const_low,d_const_low,d_const_low,d_const_low,opt_s['res_array'][0,0],opt_s['res_array'][1,0])
    #
    print(const_1_tact_full['Failure_year'],'FY 0.5',const_1_tact_med['Failure_year'],'FY 0.35',const_1_tact_low['Failure_year'],'FY 0.2')
# #----------------------------------------------------------------------------------------------
    if plots: # needs opt sim, mosaic needs pp1
        
        if present_plots:
            fig = plt.figure(figsize=(10,7))
            ax = fig.add_subplot(211)

            low_stop = 22
            mid_stop = 24
            full_stop = 26
            opt_stop = 28
            
            ax.axvline(low_stop-1,color='g',linestyle='--')
            ax.axvline(mid_stop-1,color='y',linestyle='--')
            ax.axvline(full_stop-1,color='r',linestyle='--')
            ax.axvline(opt_stop-1,color='k',linestyle='--')

            ax.plot(range(len(d_const_low)),2*d_const_low,label='Fungicide 1, 2; low dose',color='g')
            ax.plot(range(len(d_const_med)),2*d_const_med,label='Fungicide 1, 2; medium dose',color='y')
            ax.plot(range(len(d_const_full)),2*d_const_full,label='Fungicide 1, 2; full dose',color='r')
            ax.plot(range(dose_opt_int.shape[1]),2*dose_opt_int[0,:],label='Fungicide 1, 2; optimal',color='k')
            
            ax.scatter(range(low_stop),2*d_const_low[0:low_stop],color='g') # or full length eg ax.scatter(range(len(d_const_low)),2*d_const_low,color='g')
            ax.scatter(range(mid_stop),2*d_const_med[0:mid_stop],color='y')
            ax.scatter(range(full_stop),2*d_const_full[0:full_stop],color='r')
            ax.scatter(range(opt_stop),2*dose_opt_int[0,0:opt_stop],color='k')

            ax.set_ylabel('Fraction of Max. Legal Dose',fontsize=12)
            ax.set_ylim((0,1.1))
            ax.set_xlim((0,len(d_const_low)))
            ax.grid()
            
            ax2 = fig.add_subplot(212)
            opt, full, med, low = [np.zeros(len(d_const_low)) for i in range(4)]

            opt[0]  = int_1_tact['Yield_vec'][0]
            full[0] = const_1_tact_full['Yield_vec'][0]
            med[0]  = const_1_tact_med['Yield_vec'][0]
            low[0]  = const_1_tact_low['Yield_vec'][0]
            for ii in range(len(d_const_low)-1):
                if int_1_tact['Yield_vec'][ii+1]>95:
                    opt[ii+1]  = opt[ii] + int_1_tact['Yield_vec'][ii+1]
                else:
                    opt[ii+1]  = opt[ii]
                if const_1_tact_full['Yield_vec'][ii+1]>95:           
                    full[ii+1] = full[ii] + const_1_tact_full['Yield_vec'][ii+1]
                else:
                    full[ii+1] = full[ii]
                if const_1_tact_med['Yield_vec'][ii+1]>95:
                    med[ii+1]  = med[ii]  + const_1_tact_med['Yield_vec'][ii+1]
                else:
                    med[ii+1]  = med[ii]
                if const_1_tact_low['Yield_vec'][ii+1]>95:            
                    low[ii+1]  = low[ii]  + const_1_tact_low['Yield_vec'][ii+1]
                else:
                    low[ii+1]  = low[ii]
            
            ax2.axvline(low_stop-1,color='g',linestyle='--')
            ax2.axvline(mid_stop-1,color='y',linestyle='--')
            ax2.axvline(full_stop-1,color='r',linestyle='--')
            ax2.axvline(opt_stop-1,color='k',linestyle='--')

            ax2.plot(range(len(const_1_tact_low['Yield_vec'])),low,color='g')
            ax2.plot(range(len(const_1_tact_med['Yield_vec'])),med,color='y')
            ax2.plot(range(len(const_1_tact_full['Yield_vec'])),full,color='r')
            ax2.plot(range(len(int_1_tact['Yield_vec'])),opt,color='k')
            
            ax2.scatter(range(low_stop),low[0:low_stop],color='g') # or full length e.g. # ax2.scatter(range(len(const_1_tact_low['Yield_vec'])),low,color='g')
            ax2.scatter(range(mid_stop),med[0:mid_stop],color='y')
            ax2.scatter(range(full_stop),full[0:full_stop],color='r')
            ax2.scatter(range(opt_stop),opt[0:opt_stop],color='k')

            ax2.set_ylabel('Yield',fontsize=12)
            ax2.grid()

            L,R,B,T = (0.12,0.9,0.12,0.9)
            fig.tight_layout(rect=[0, 0.03, 1, 0.95])
            fig.subplots_adjust(left=L,right=R,bottom=B,top=T)
            xlab_ax = fig.add_axes([L, 0.08, R-L-0.06, 0.01])
            def inv_ax(big_ax):
                for jj in ['top','bottom','right','left']:
                    big_ax.spines[jj].set_color('none')
                big_ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
                return big_ax
            xlab_ax = inv_ax(xlab_ax)
            xlab_ax.set_xlabel('Season Number',fontsize=16)
            ax2.set_xlim((0,len(d_const_low)))

            fig.legend()    
# #----------------------------------------------------------------------------------------------
        if season_plt:
            res_freq_ratio = opt_s['res_array'][0,:]
            res_freq_ratio = [opt_s['res_array'][0,i]/(opt_s['res_array'][1,i]+opt_s['res_array'][0,i]) for i in range(len(res_freq_ratio))]
            res_freq_ratio = res_freq_ratio[0:(opt_s['success_years']+1)]
            print(res_freq_ratio)
            Season_plotter('r',SRV1=not_int_1_tact['Yield_vec'],SRV2=int_1_tact['Yield_vec'],RV1=res_freq_ratio,RV1_y_lim=(0,1),SRV1_label='Y1',SRV2_label='Y1_int',RV1_label='RF_ratio')
# #----------------------------------------------------------------------------------------------
        for i in range(0,opt_s['success_years'],2): # [0,2,3]:#
            if overlay_plt:
                if interp:
                    Overlay_plotter(Col_mat= opt_s['H_int'][:,:,i],Col_label='G',Col_bds_vec=(1,),Con_mats=(opt_s['Y_int'][:,:,0,i]),Con_levels=([params.Yield_threshold]),Con_inline=('inline'),title='G in year %s' % i)
                else:
                    Overlay_plotter(Col_mat= opt_s['H'][:,:,i],Col_label='G',Col_bds_vec=(1,),Con_mats=(opt_s['Yield_array'][:,:,0,i]),Con_levels=([params.Yield_threshold]),Con_inline=('inline'),title='G in year %s' % i)
            if contour_plt:
                if g_or_y == 'G':
                    Overlay_plotter(Con_mats=(opt_s['Yield_array'][:,:,0,i],opt_s['H'][:,:,i]),Con_levels=([params.Yield_threshold],[1,3,4,4.5,5,5.5]),Con_inline=('inline','inline'),title='G in year %s' % i)
                else:
                    if n_d_i==n_d:
                        contour = np.linspace(np.amin(opt_s['Y_FD_by_dose_to_plot'][:,:,i]),np.amax(opt_s['Y_FD_by_dose_to_plot'][:,:,i]),41)
                        Overlay_plotter(Con_mats=(opt_s['Yield_array'][:,:,0,i],opt_s['Y_FD_by_dose_to_plot'][:,:,i]),Con_levels=([params.Yield_threshold],contour),Con_inline=('inline','inline'),title='G in year %s' % i)
                        # print(Y_FD_by_dose)
# #----------------------------------------------------------------------------------------------
        if surface_plt:
            X, Y = np.meshgrid(np.linspace(0,1,n_d),np.linspace(0,1,n_d))
            for i in [0,2,4]:
                fig2 = plt.figure()
                ax2 = fig2.gca(projection='3d')
                surf2 = ax2.plot_wireframe(X,Y,opt_s['H'][:,:,i])
                ax2.set_title('G')
# #----------------------------------------------------------------------------------------------
        if mosaic_plt:
            if g_or_y == 'G':
                optimal_mosaic(opt_s['H_int'],opt_s['Y_int'],dose_opt_int,opt_s['success_years'],last_9=True,type_plot='contour',scatter_on=False)
                # optimal_mosaic(opt_s['H_int'],opt_s['Y_int'],dose_opt_int,opt_s['success_years'],last_9=True,type_plot='contourf')
            if g_or_y == 'Y':
                ZZ = np.zeros(opt_s['Y_FD_by_dose_to_plot'].shape)
                optimal_mosaic(opt_s['Y_FD_by_dose_to_plot'],opt_s['Y_int'],dose_opt_int,opt_s['success_years'],last_9=True,type_plot='contour',bounds_vector='linspace')
                optimal_mosaic(ZZ,opt_s['Y_int'],dose_opt_int,opt_s['success_years'],last_9=True,type_plot='contour',bounds_vector='linspace')
                optimal_mosaic(opt_s['Y_FD_by_dose_to_plot'],opt_s['Y_int'],dose_opt_int,opt_s['success_years'],last_9=True,type_plot='contourf',bounds_vector='linspace')
            if g_or_y == 'D':
                ZZ = np.zeros(pp_1['Diff'].shape)
                optimal_mosaic(pp_1['Diff'],opt_s['Y_int'],dose_opt_int,opt_s['success_years'],last_9=False,type_plot='contour',bounds_vector='linspace')
                optimal_mosaic(ZZ,opt_s['Y_int'],dose_opt_int,opt_s['success_years'],last_9=False,type_plot='contour',bounds_vector='linspace')
                optimal_mosaic(pp_1['Diff'],opt_s['Y_int'],dose_opt_int,opt_s['success_years'],last_9=False,type_plot='contourf',bounds_vector='linspace')
                for i in range(opt_s['Y_int'].shape[-1]): # n years
                    Con_attrs2 = {'Con_mats': (opt_s['Y_int'][:,:,0,i]),'Con_levels':([0,95]),'contourf':0}#'k'
                    Overlay_plotter(Con_attr=Con_attrs2,figure_attributes={'title':'Year %s'%i})
                    print(opt_s['Y_int'][:,:,0,i])
# #----------------------------------------------------------------------------------------------
## end of 'if plots:'
# #----------------------------------------------------------------------------------------------
    if phi_cont_plt:  # needs opt sim
        n = 8
        dose = 0.5
        C = 0.1 # 0.13
        ##
        YieldV_p_c_plt = np.zeros((n,n))
        xy_array  = np.zeros((n,n))
        phi_1 = np.linspace(0,1,n)
        phi_2 = np.linspace(0,1,n)
        ##
        for i in range(len(phi_1)):
            for j in range(len(phi_2)):
                p_c_plt_output = master_loop_one_tactic([dose],[dose],[dose],[dose],phi_1[i],phi_2[j])
                YieldV_p_c_plt[i,j] = p_c_plt_output['Yield_vec']
                xy_array[i,j] = C - phi_1[i]*phi_2[j]
        ##
        Rv1=int_1_tact['Res_vec_1'][0:-1]
        Rv2=int_1_tact['Res_vec_2'][0:-1]
        
        markers = ['+']*n_seasons
        for i in range(max(0.5*len(markers),len(markers)-5),len(markers)):
            markers[i] = '$' + str(i) + '$'
        
        cont = [params.Yield_threshold,99,99.05,99.06,99.0668,99.07]
        cont = [0,params.Yield_threshold]
        Con_attrs = {'Con_mats': (YieldV_p_c_plt),'Con_levels':(cont),'Con_inline':(None),'Contourf':('grey')}
        Scat_att  = {'x_scat':Rv1,'y_scat':Rv2,'markers':markers}
        Fig_att = {'x_lab':r'Resistance Frequency 1 ($\phi_1$)','y_lab':r'Resistance Frequency 2 ($\phi_2$)','grid':False}
        Overlay_plotter(Con_attr=Con_attrs,Scat_attr=Scat_att,figure_attributes=Fig_att)

        # if g_or_y == 'G':
        #     Overlay_plotter(Con_mats=(YieldV_p_c_plt,xy_array),Con_inline=('inline',None),Con_levels=(cont,[0]),x_scat=Rv1,y_scat=Rv2,title=r'Yield, $\phi_1$, $\phi_2$',x_lab=r'$\phi_1$',y_lab=r'$\phi_2$')
        #     Overlay_plotter(Con_mats=(YieldV_p_c_plt,xy_array),Con_inline=('inline',None),Con_levels=(cont,[0]),title=r'Yield, $\phi_1$, $\phi_2$',x_lab=r'$\phi_1$',y_lab=r'$\phi_2$')
# #----------------------------------------------------------------------------------------------
    if c_cont_plt: # needs opt sim
        # C_vec = np.concatenate(([10**(-8),10**(-6)],np.linspace(10**(-4),0.2,n_c-2)))

        # a = [10**(-6),10**(-3),4*10**(-3),8*10**(-3),1.2*10**(-2),1.6*10**(-2)]
        # b = np.linspace(2*10**(-2),0.5,n_p-6)
        # p = np.concatenate((a,b))

        phi_v = np.zeros((2,len(int_1_tact['Res_vec_1'])))
        phi_v[0,:] = int_1_tact['Res_vec_1']
        phi_v[1,:] = int_1_tact['Res_vec_2']
        phi_v=phi_v[:,0:-1] # ignore the final point
        c_bottom = phi_v[0,0]*phi_v[1,0]
        c_top = phi_v[0,-1]*phi_v[1,-1]
        c_top = c_top + (1-c_top)*0.1
        print(phi_v,'phi_v')
        C_contours(c_bottom=c_bottom,c_top=c_top,phi_vec=phi_v)
# #----------------------------------------------------------------------------------------------
    # default is: animate_not_int = 0, animate_int = 1
    if animate:
        optimal_animator(not_int_1_tact['Res_vec_1'],not_int_1_tact['Res_vec_2'],int_1_tact['Res_vec_1'],int_1_tact['Res_vec_2'],not_int_1_tact['Yield_vec'],int_1_tact['Yield_vec'],opt_s['Yield_array'],opt_s['Y_int'],opt_s['H'],opt_s['H_int'],dose_opt_not_int,dose_opt_int,interp,opt_s['success_years'],n_seasons)
# #----------------------------------------------------------------------------------------------
plt.show()