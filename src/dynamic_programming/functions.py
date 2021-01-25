import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from math import ceil, log, floor, log10, exp
from scipy.optimize import fsolve
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import RegularGridInterpolator
import warnings

from utils.plotter_HRHR    import Overlay_plotter, colourmap_fn
from utils.functions_HRHR  import interpolate, master_loop_one_tactic, master_loop_grid_of_tactics
from utils.parameters_HRHR import PARAMS


# #----------------------------------------------------------------------------------------------        
def Int_Y_FD(n_phi = 2,dose = 0.5,output_size=10):
    Yield_vec = np.zeros((n_phi,n_phi))
    phi_1, phi_2 = [np.linspace(0,1,n_phi) for i in range(2)]
    # phi_2 = np.linspace(0,1,n_phi)
    for i in range(n_phi):
        for j in range(n_phi):
            output = master_loop_one_tactic([dose],[dose],[dose],[dose],phi_1[i],phi_2[j]) # but full dose might be 0.5
            Yield_vec[i,j] = output['Yield_vec']
    Interpolate_Y_at_FD, Int_Y_at_FD = interpolate(Yield_vec,len(Yield_vec[:,0]),output_size,'linear')
    return Interpolate_Y_at_FD, Int_Y_at_FD







# #----------------------------------------------------------------------------------------------        
# currently only correct for sexual reproduction
def optimal_simulator(r1, r2, n_doses, output_size, n_seasons, interp, 
        Interpolate_Y_at_FD=None, always_on_contour=0, interp_type=None, 
        equal_dose_on_contour=0, G_or_yield=None, n_phi=None, Diff=None):
    
    # default settings:
    n_seas = 1
    upper_Y = 100
    upper_G = 100
    upper_P = 100

    ##
    if interp_type is None:
        interp_type = 'linear'
    if G_or_yield is None:
        G_or_yield = 'Y'
    if n_phi is None:
        n_phi = n_doses
    # #----------------------------------------------------------------------------------------------        
    # initialise
    res_array                 = 2*np.ones((2,n_seasons+1))
    dose_array,dose_array_int = [2*np.ones((2,n_seasons)) for i in range(2)]
    ##
    Yield_array              = np.zeros((n_doses,n_doses,n_seas,n_seasons))
    S1,S2 = [np.zeros((n_doses,n_doses,n_seasons+1)) for i in range(2)]
    ##
    G,H = [upper_G*np.ones((n_doses,n_doses,n_seasons+1)) for i in range(2)]
    ##
    H_int,G_int,G_int_lin = [upper_G*np.ones((output_size,output_size,n_seasons)) for i in range(3)] #+1
    ##
    Y_int                    = np.zeros((output_size,output_size,n_seas,n_seasons))
    Y_FD_by_dose,Y_FD_by_dose_complete = [np.zeros((n_phi,n_phi,n_seasons)) for i in range(2)]
    ##
    Y_FD_by_dose_to_plot     = np.zeros((output_size,output_size,n_seasons))
    P_speed_by_dose,P_speed_by_dose_complete = [upper_P*np.ones((n_phi,n_phi,n_seasons)) for i in range(2)]
    ##
    res_array[0,0] = r1
    res_array[1,0] = r2
    # #----------------------------------------------------------------------------------------------        
    # Yield at full dose as a function of phi
    if G_or_yield == 'Y':
        if Interpolate_Y_at_FD is None:
            Interpolate_Y_at_FD = Int_Y_FD(n_phi = n_phi)[0] # , Int_Y_at_FD
    # #----------------------------------------------------------------------------------------------        
    # speed in phi-space
    if G_or_yield== 'D':
        Interpolate_Phi_Speed = interpolate(Diff,len(Diff[:,0]),output_size,'linear')[0]
        # Interpolate_Phi_Speed, Int_Phi_Speed = interpolate(Diff,len(Diff[:,0]),output_size,'linear')
    # #----------------------------------------------------------------------------------------------        
    for k in range(n_seasons):
        if k >= 1:
            if Yield_array[n_doses-1,n_doses-1,0,k-1] > PARAMS.Yield_threshold: # survived last year
                grid_output = master_loop_grid_of_tactics(n_doses,n_seas,res_array[0,k],res_array[1,k])#y_t
                Yield_array[:,:,:,k] = grid_output['Yield']
                Res_array_1          = grid_output['Res_array_1']
                Res_array_2          = grid_output['Res_array_2']
                S1                   = grid_output['Selection_array_1']
                S2                   = grid_output['Selection_array_2']
                success_years = k
        else:
            grid_output = master_loop_grid_of_tactics(n_doses,n_seas,res_array[0,k],res_array[1,k])#y_t
            Yield_array[:,:,:,k] = grid_output['Yield']
            Res_array_1          = grid_output['Res_array_1']
            Res_array_2          = grid_output['Res_array_2']
            S1                   = grid_output['Selection_array_1']
            S2                   = grid_output['Selection_array_2']
            success_years = k
        ##
        if success_years == k: # so still surviving
            if G_or_yield == 'G':
                for i in range(n_doses):
                    for j in range(n_doses):
                        H[i,j,k] = S1[i,j,1]*S2[i,j,1]
                        if Yield_array[i,j,0,k]>PARAMS.Yield_threshold:
                            G[i,j,k] = S1[i,j,1]*S2[i,j,1]
                if np.amin(G[:,:,k])<upper_G:
                    index_G = np.argwhere(G[:,:,k]==np.amin(G[:,:,k]))
                    ii = index_G[0][0]
                    jj = index_G[0][1]
                    dose_array[0,k] = 0.5*ii/(n_doses-1)
                    dose_array[1,k] = 0.5*jj/(n_doses-1)
                    res_array[0,k+1]  = Res_array_1[ii,jj,1]
                    res_array[1,k+1]  = Res_array_2[ii,jj,1]
            for i in range(n_phi):
                for j in range(n_phi):
                    if G_or_yield == 'Y':
                        Y_FD_by_dose_complete[i,j,k] = Interpolate_Y_at_FD(Res_array_1[i,j,1],Res_array_2[i,j,1])
                        if Yield_array[i,j,0,k]>PARAMS.Yield_threshold: # so tactic valid
                            Y_FD_by_dose[i,j,k] = Interpolate_Y_at_FD(Res_array_1[i,j,1],Res_array_2[i,j,1])
                        Y_FD_by_dose_int = interpolate(Y_FD_by_dose[:,:,k],n_phi,output_size,'linear')[1] # Interpolate_Y_FD_plot, Y_FD_by_dose_int, without [0] in function call
                        Y_FD_by_dose_to_plot[:,:,k] = interpolate(Y_FD_by_dose_complete[:,:,k],n_phi,output_size,'linear')[1] # Interpolate_Y_FD_plot, Y_FD_by_dose_int, without [0] in function call
                    if G_or_yield == 'D':
                        P_speed_by_dose_complete[i,j,k] = Interpolate_Phi_Speed(Res_array_1[i,j,1],Res_array_2[i,j,1])
                        if Yield_array[i,j,0,k]>PARAMS.Yield_threshold: # so tactic valid
                            P_speed_by_dose[i,j,k] = Interpolate_Phi_Speed(Res_array_1[i,j,1],Res_array_2[i,j,1])
                        P_speed_by_dose_int = interpolate(P_speed_by_dose[:,:,k],n_phi,output_size,'linear')[1] #Interpolate_P_speed_in
    # #----------------------------------------------------------------------------------------------        
            if interp==1:
                # define interpolated functions
                if G_or_yield == 'G':
                    Interpolate_H, Int_H = interpolate(H[:,:,k],H.shape[0],output_size,interp_type)
                    H_int[:,:,k] = Int_H
                    Int_G_lin = interpolate(G[:,:,k],G.shape[0],output_size,'linear')[1] #Interpolate, 
                    Int_G = interpolate(G[:,:,k],G.shape[0],output_size,interp_type)[1] #Interpolate, 
                    G_int[:,:,k] = Int_G
                    G_int_lin[:,:,k] = Int_G_lin
                Interpolate_Y, Int_Y = interpolate(Yield_array[:,:,0,k],Yield_array.shape[0],output_size,interp_type)
                Y_int[:,:,0,k] = Int_Y
                Interpolate_R1 = interpolate(Res_array_1[:,:,1],Res_array_1.shape[0],output_size,interp_type)[0] # , Int_R1
                Interpolate_R2 = interpolate(Res_array_2[:,:,1],Res_array_2.shape[0],output_size,interp_type)[0] # , Int_R2
                # #----------------------------------------------------------------------------------------------                    
                # full dose or straight contours # was designed in case max wasn't in top corner but was for high doses
                if G_or_yield == 'G':
                    matrix = G_int_lin[:,:,k]
                    Metric_1_FD = np.amin(np.concatenate((G_int_lin[:,0,k],G_int_lin[0,:,k],G_int_lin[:,-1,k],G_int_lin[-1,:,k])))
                    Metric_2_FD = np.amin(G_int_lin[:,:,k])
                if G_or_yield == 'Y':
                    matrix = Y_FD_by_dose_int
                    Metric_1_FD = np.amax(np.concatenate((Y_FD_by_dose_int[:,0],Y_FD_by_dose_int[0,:],Y_FD_by_dose_int[:,-1],Y_FD_by_dose_int[-1,:])))
                    Metric_2_FD = np.amax(Y_FD_by_dose_int)
                if G_or_yield=='D':
                    matrix = P_speed_by_dose_int
                    Metric_1_FD = np.amin(np.concatenate((P_speed_by_dose_int[:,0],P_speed_by_dose_int[0,:],P_speed_by_dose_int[:,-1],P_speed_by_dose_int[-1,:])))
                    Metric_2_FD = np.amin(P_speed_by_dose_int)
                # #----------------------------------------------------------------------------------------------                    
                index = np.argwhere(matrix==Metric_1_FD)
                i_y_chosen, j_y_chosen = index[0] # corresponds to the pair found
                # #----------------------------------------------------------------------------------------------                    
                if always_on_contour == 0 and Metric_1_FD==Metric_2_FD and success_years == k: # so still surviving
                    # potentially should switch a tiny bit earlier if there exists a minimum elsewhere on 95% contour
                    # use linear interpolation for this condition, even if using cubic elsewhere.
                    dose_array_int[0,k] = 0.5*i_y_chosen/(output_size-1)
                    dose_array_int[1,k] = 0.5*j_y_chosen/(output_size-1)
                    ##
                    res_array[0,k+1] = Interpolate_R1(2*dose_array_int[0,k],2*dose_array_int[1,k])
                    res_array[1,k+1] = Interpolate_R2(2*dose_array_int[0,k],2*dose_array_int[1,k])
                    print(res_array[0,k+1],res_array[1,k+1],k,'str',G_or_yield)
                # #----------------------------------------------------------------------------------------------        
                # here downwards selects point on the contour: x_val_chosen, y_val_chosen
                else:
                    def equations(y):
                        return Interpolate_Y(x_val,y)-PARAMS.Yield_threshold # x_val is in [0,1]
                    ##
                    cont_min_exact = upper_Y
                    cont_max_exact = 0
                    ig = 0.2
                    x_val_chosen, y_val_chosen = (2,2) # in case nothing on contour works
                    ##
                    if equal_dose_on_contour==1:
                        def equations_equal(y):
                            return Interpolate_Y(y,y)-PARAMS.Yield_threshold # x_val is in [0,1]
                        y_val = fsolve(equations_equal,ig)
                        x_val_chosen = y_val
                        y_val_chosen = y_val
                    ##
                    if equal_dose_on_contour==0:
                        for i in range(output_size):
                            x_val = i/(output_size-1) # in [0,1]
                            y_val = fsolve(equations,ig) # in [0,1]          # given x, find y such that are on 95% contour
                            if y_val == ig:
                                y_val = fsolve(equations,0.5)
                            
                            # print(x_val,y_val)
                            if y_val<1 and y_val>0: # <= >=
                                if G_or_yield == 'G':
                                    if Interpolate_H(x_val,y_val) < cont_min_exact:
                                        cont_min_exact = Interpolate_H(x_val,y_val)
                                        x_val_chosen , y_val_chosen = x_val, y_val # on the contour, in [0,1]
                                if G_or_yield == 'Y': ####
                                    p1 = Interpolate_R1(x_val,y_val)
                                    p2 = Interpolate_R2(x_val,y_val)
                                    if Interpolate_Y_at_FD(p1,p2) > cont_max_exact:
                                        cont_max_exact = Interpolate_Y_at_FD(p1,p2)
                                        x_val_chosen , y_val_chosen = x_val, y_val # on the contour, in [0,1]
                                if G_or_yield == 'D':
                                    p1 = Interpolate_R1(x_val,y_val)
                                    p2 = Interpolate_R2(x_val,y_val)
                                    if Interpolate_Phi_Speed(p1,p2) < cont_min_exact:
                                        cont_min_exact = Interpolate_Phi_Speed(p1,p2)
                                        x_val_chosen , y_val_chosen = x_val, y_val # on the contour, in [0,1]
                    
                    
                    # #----------------------------------------------------------------------------------------------                            
                    dose_array_int[0,k] = 0.5*x_val_chosen        # in [0,0.5]
                    dose_array_int[1,k] = 0.5*y_val_chosen #0.5*contour_dose_interp[i_chosen,k] #  in [0,0.5], corresponds to the contour
                    
                    res_array[0,k+1] = Interpolate_R1(x_val_chosen,y_val_chosen) # flip?
                    res_array[1,k+1] = Interpolate_R2(x_val_chosen,y_val_chosen) # flip?
                    print('Doses',dose_array_int[0,k],dose_array_int[1,k],dose_array_int[0,k]/(dose_array_int[0,k]+dose_array_int[1,k]),'D ratio RF',res_array[0,k+1]/(res_array[0,k+1]+res_array[1,k+1]),k,'CD',G_or_yield,res_array[0,k+1],res_array[1,k+1],'RFs')
    
    dictionary = {'H': H,
        'G': G,
        'Yield_array': Yield_array,
        'dose_array': dose_array,
        'res_array': res_array,
        'success_years': success_years,
        'H_int': H_int,
        'G_int': G_int,
        'Y_int': Y_int,
        'dose_array_int': dose_array_int,
        'Y_FD_by_dose': Y_FD_by_dose,
        'Y_FD_by_dose_to_plot': Y_FD_by_dose_to_plot
        }
    
    return dictionary










# #----------------------------------------------------------------------------------------------      
def optimal_mosaic(H_int,Y_int,dose,success_years,last_9=False,type_plot=None,bounds_vector=None,scatter_on=True):
    

    ##
    label = 'G'
    if bounds_vector is None:
        bounds_vector = (1,)
    
    numberofdoses_H = H_int.shape[0]
    numberofdoses_Y = Y_int.shape[0]
    ##
    x_H = np.linspace(0-0.5/(numberofdoses_H-1),1+0.5/(numberofdoses_H-1),numberofdoses_H+1)
    y_H = np.linspace(0-0.5/(numberofdoses_H-1),1+0.5/(numberofdoses_H-1),numberofdoses_H+1)
    X_H, Y_H = np.meshgrid(x_H,y_H)
    ##
    # x_Y = np.linspace(0-0.5/(numberofdoses_Y-1),1+0.5/(numberofdoses_Y-1),numberofdoses_Y+1)
    # y_Y = np.linspace(0-0.5/(numberofdoses_Y-1),1+0.5/(numberofdoses_Y-1),numberofdoses_Y+1)
    # X_Y, Y_Y = np.meshgrid(x_Y,y_Y)
    ##
    labels = ([0,PARAMS.Yield_threshold])
    alpha = 1
    ##
    x1_H = np.linspace(0,1,numberofdoses_H)
    y1_H = np.linspace(0,1,numberofdoses_H)
    X1_H, Y1_H = np.meshgrid(x1_H,y1_H)
    ##
    x1_Y = np.linspace(0,1,numberofdoses_Y)
    y1_Y = np.linspace(0,1,numberofdoses_Y)
    X1_Y, Y1_Y = np.meshgrid(x1_Y,y1_Y)
    ##
    top = min(success_years+1,len(Y_int[0,0,0,:]))
    jump = ceil((top+1)/9)
    div = ceil(top/jump)
    x = ceil(div**(0.5))
    fig, shared_ax = plt.subplots(x, x, sharex=True, sharey=True, figsize=(10,9))

    j = 0
    if last_9:
        bottom = max(top-9,0)
        vector = range(bottom,top,1)
    else:
        vector = range(0,top,jump)
    for i in vector:
        # sub = str(x) + str(x) + str(j)
        ii = j % x
        jj = floor(j/x)
        
        j = j+1
        ax = shared_ax[jj,ii]

        print(np.ndim(H_int))
        if np.ndim(H_int) == 3:
            matrix = H_int[:,:,i]
        else:
            matrix = H_int
        array = Y_int[:,:,0,i]

        title = 'Year %s' % i

        if bounds_vector is None:
            Maximum = ceil(np.amax(matrix))
            bounds    = np.linspace(bounds_vector[0],Maximum,Maximum-bounds_vector[0]+1)
            ticks    = np.linspace(bounds_vector[0],Maximum,Maximum-bounds_vector[0]+1)
        else:
            Minimum = floor(np.amin(matrix))
            Maximum = ceil(np.amax(matrix))
            bounds = np.linspace(Minimum,Maximum,41)#(Maximum-Minimum)*8+1
            ticks    = np.linspace(Minimum,Maximum,Maximum+1-Minimum)
        
        cmap, norm = colourmap_fn(bounds)

        if type_plot is None:
            im = ax.pcolormesh(X_H,Y_H,np.transpose(matrix),cmap=cmap, norm=norm,alpha = alpha,linewidth=0,antialiased=True)# np.transpose(matrix)
        if type_plot == 'contour':
            im = ax.contour(X1_H,Y1_H,np.transpose(matrix),levels=bounds)
            ax.clabel(im,inline=1,fontsize=9)
        if type_plot == 'contourf':
            im = ax.contourf(X1_H,Y1_H,np.transpose(matrix),cmap=cmap,levels=bounds)
        # antialiased - attempt to remove grid lines that appear due to overlap of alpha=0.5 squares

        ax.contourf(X1_Y,Y1_Y,np.transpose(array),colors=('k'),alpha = 0.5,levels = labels)# np.transpose(matrix)
        ax.set_xlim((0,1))
        ax.set_ylim((0,1))
        
        if scatter_on:
            if i <success_years:
                ax.scatter(2*dose[0,i],2*dose[1,i],color='c',clip_on=False)

        if type_plot !='contour':
            divider = make_axes_locatable(ax)
            cax = divider.new_horizontal(size="5%", pad=0.3, pack_start=False)
            # create a second axes for the colorbar
            ax2 = fig.add_axes(cax)
            mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=ticks, boundaries=bounds)#, format='%1i')
            ax2.set_ylabel(label)

    ##
        ax.set_title(title,fontsize=18)
        if jj == 2:
            ax.set_xlabel(r'$F_1$ dose',fontsize=16)
            ax.set_xticks([0,0.2,0.4,0.6,0.8,1])
        if ii == 0:
            ax.set_ylabel(r'$F_2$ dose',fontsize=16)
            ax.set_yticks([0,0.2,0.4,0.6,0.8,1])

    mpl.pyplot.subplots_adjust(wspace=0.5,hspace=0.5)
    return None








# #----------------------------------------------------------------------------------------------      
def phi_plane_jumps(C_cont=None,n_steps = 3,n_phi = 12,n_d=10,output_size = 21,phi_grid=None):
    dose = 0.5
    ig = 0.2
    n_s = 2
    ##
    start_val = 2
    start_val_sr = 0.9
    start_val_Y = 50
    ##
    Yield_vec    = np.zeros((n_phi,n_phi))
    phi_1, phi_2 = [np.linspace(0,1,n_phi) for i in range(2)]
    SR1, SR2     = [start_val_sr*np.ones((n_d,n_d,n_steps)) for i in range(2)]
    dist1, dist2 = [start_val*np.ones((n_d,n_d,n_steps))    for i in range(2)]
    Rv1, Rv2     = [start_val*np.ones((n_d,n_d,n_steps))    for i in range(2)]
    ##
    Diff  = np.zeros((n_phi,n_phi))
    phi_1_values,phi_2_values = [np.zeros(n_steps) for i in range(2)]
    ##
    # generate phi plane for full dose
    for i in range(len(phi_1)):
        for j in range(len(phi_2)):
            # Res_vec_1, Res_vec_2, PRR, PRS, PSR, PSS, Yield_vec[i,j], Innoc_vec,  Selection_array1, Selection_array2, Failure_year, Sol_array, t_vec = master_loop_one_tactic([dose],[dose],[dose],[dose],phi_1[i],phi_2[j])
            one_t_output = master_loop_one_tactic([dose],[dose],[dose],[dose],phi_1[i],phi_2[j])
            Yield_vec[i,j] = one_t_output['Yield_vec']
    #
    Interpolate_Y_FD = interpolate(Yield_vec,n_phi,output_size)[0] # , Int_Y_FD
    # #----------------------------------------------------------------------------------------------      
    if phi_grid is None:
        Yield = start_val_Y*np.ones((n_d,n_d,n_s,n_steps))
        Y2    = start_val_Y*np.ones((n_d,n_d,n_steps))
        Res_array_1, Res_array_2 = [np.zeros((n_d,n_d,n_s+1,n_steps)) for i in range(2)]
        def equations_FD_phi_2_min(y):
            return Interpolate_Y_FD(1,y) - C_cont
        def equations_FD_equal(y):
            return Interpolate_Y_FD(y,y) - C_cont
        def equations_FD_phi_1(x):
            return Interpolate_Y_FD(x,y_val) - C_cont
    ###
        phi_2_min = fsolve(equations_FD_phi_2_min,ig)
        if phi_2_min > 1 or phi_2_min < 0:
            warnings.warn('Warning! Min out of bounds, ',phi_2_min,' corrected to 10**(-4)')
            phi_2_min = 10**(-4)

        phi_2_max = fsolve(equations_FD_equal,ig)
        if phi_2_max > 1 or phi_2_max < 0:
            phi_2_max = 1
            warnings.warn('Warning! Max out of bounds')
        ##
        phi_2_values = np.linspace(phi_2_min,phi_2_max,n_steps) # gives us n_steps values along contour in lower triangle
        for i in range(n_steps):
            y_val = phi_2_values[i]
            phi_1_values[i] = min(fsolve(equations_FD_phi_1,ig),1)        
    ###
        for k in range(n_steps):
            # LTY, TY, FY, Yield[:,:,:,k], Res_array_1[:,:,:,k], Res_array_2[:,:,:,k], PRR_array, PRS_array, PSR_array, PSS_array, Selection_array_1, Selection_array_2, Innoc_array = master_loop_grid_of_tactics(n_d,n_s,phi_1_values[k],phi_2_values[k],yield_stopper=PARAMS.Yield_threshold) # doses
            grid_output = master_loop_grid_of_tactics(n_d,n_s,phi_1_values[k],phi_2_values[k],yield_stopper=PARAMS.Yield_threshold) # doses
            Yield[:,:,:,k]       = grid_output['Yield']
            Res_array_1[:,:,:,k] = grid_output['Res_array_1']
            Res_array_2[:,:,:,k] = grid_output['Res_array_2']


            for i in range(n_d):
                for j in range(n_d):
                    if Yield[i,j,0,k]>PARAMS.Yield_threshold: # so acceptable dose
                        Rv1[i,j,k]   = Res_array_1[i,j,1,k]
                        Rv2[i,j,k]   = Res_array_2[i,j,1,k]
                        dist1[i,j,k] = Res_array_1[i,j,1,k]-phi_1_values[k]
                        dist2[i,j,k] = Res_array_2[i,j,1,k]-phi_2_values[k]
                        Y2[i,j,k] = Interpolate_Y_FD(Res_array_1[i,j,1,k],Res_array_2[i,j,1,k])
                        SR1[i,j,k] = Rv1[i,j,k]/phi_1_values[k]
                        SR2[i,j,k] = Rv2[i,j,k]/phi_2_values[k]
    # #----------------------------------------------------------------------------------------------      
    else:
        Yield = start_val_Y*np.ones((n_d,n_d,n_s,n_phi,n_phi))
        Y2    = start_val_Y*np.ones((n_d,n_d,n_phi,n_phi))
        Res_array_1 = np.zeros((n_d,n_d,n_s+1,n_phi,n_phi))
        Res_array_2 = np.zeros((n_d,n_d,n_s+1,n_phi,n_phi))
        phi_1_values = np.linspace(0,1,n_phi)
        phi_2_values = np.linspace(0,1,n_phi)
        for ii in range(n_phi):
            for jj in range(n_phi):
                if Interpolate_Y_FD(phi_1[ii],phi_2[jj])>PARAMS.Yield_threshold: # and we use yield threshold on the grid
                    # LTY, TY, FY, Yield[:,:,:,ii,jj], Res_array_1[:,:,:,ii,jj], Res_array_2[:,:,:,ii,jj], PRR_array, PRS_array, PSR_array, PSS_array, Selection_array_1, Selection_array_2, Innoc_array = master_loop_grid_of_tactics(n_d,n_s,phi_1[ii],phi_2[jj],yield_stopper=PARAMS.Yield_threshold) # doses
                    grid_2_output = master_loop_grid_of_tactics(n_d,n_s,phi_1[ii],phi_2[jj],yield_stopper=PARAMS.Yield_threshold) # doses
                    Yield[:,:,:,ii,jj]       = grid_2_output['Yield']
                    Res_array_1[:,:,:,ii,jj] = grid_2_output['Res_array_1']
                    Res_array_2[:,:,:,ii,jj] = grid_2_output['Res_array_2']

                    for i in range(n_d):
                        for j in range(n_d):
                            if Yield[i,j,0,ii,jj]>PARAMS.Yield_threshold: # so acceptable dose
                                Y2[i,j,ii,jj] = Interpolate_Y_FD(Res_array_1[i,j,1,ii,jj],Res_array_2[i,j,1,ii,jj])
                Diff[ii,jj] = Interpolate_Y_FD(phi_1[ii],phi_2[jj]) - np.amax(Y2[:,:,ii,jj])
    # #----------------------------------------------------------------------------------------------      
    dictionary = {'Rv1': Rv1, 'Rv2': Rv2, 'Res_array_1': Res_array_1, 'Res_array_2': Res_array_2, 'dist1': dist1, 'dist2': dist2, 'Yield': Yield, 'Y2': Y2, 'SR1': SR1, 'SR2': SR2, 'phi_1_values': phi_1_values, 'phi_2_values': phi_2_values, 'Yield_vec': Yield_vec, 'Diff': Diff, 'Interpolate_Y_FD': Interpolate_Y_FD}
    return dictionary




# #----------------------------------------------------------------------------------------------      
def transformed_yield_contours(grid=False,n_phi=5,n_steps=5,dose_descriptors=False,Interpolate_Y_FD=None,phi_vector=None,Y_cont=None,n_doses=5,x_lim=None,x_points=None,y_points=None,colour_vec=None,n_s=None,hlines=None,bin_from=None,log_ind = True,title=True):
    fig = plt.figure(figsize=(12,6))
    ax  = fig.add_subplot(111)
    cmap = plt.get_cmap('jet')
    col = 1

    
    if phi_vector is None:
        phi_vector = range(n_steps)
    if Interpolate_Y_FD is not None:
        Interp_Y_FD = Interpolate_Y_FD
        No_R_pt = Interpolate_Y_FD(0,0)
        if log_ind:
            ax.fill_between([log10(No_R_pt - PARAMS.Yield_threshold),10],-10,10, facecolor='grey')
            # ax.axvline(No_R_pt,linestyle='--',color='k') # resistance free yield
            # ax.axvline(log10(No_R_pt - PARAMS.Yield_threshold),linestyle='--',color='k')
    ##
    
    # ax.axvline(No_R_pt,linestyle='--',color='k') # resistance free yield
    if not log_ind:
        ax.axvline(PARAMS.Yield_threshold,linestyle='--',color='k')


    ##
    if x_points is None and y_points is None:
        if grid:
            # Rv1, Rv2, Res_A1, Res_A2, dist1, dist2, Yield, Y2, SR1, SR2, phi_1_values, phi_2_values, Yield_vec, Diff, Interpolate_Y_FD = phi_plane_jumps(phi_grid=1,n_phi=n_phi,n_d=n_doses,n_steps=n_steps)
            pp_jumps_out = phi_plane_jumps(phi_grid=1,n_phi=n_phi,n_d=n_doses,n_steps=n_steps)
            Y2           = pp_jumps_out['Y2']
            Res_A1       = pp_jumps_out['Res_array_1']
            Res_A2       = pp_jumps_out['Res_array_2']
            phi_1_values = pp_jumps_out['phi_1_values']
            phi_2_values = pp_jumps_out['phi_2_values']
            if Interpolate_Y_FD is None:
                Interp_Y_FD = pp_jumps_out['Interpolate_Y_FD']
                No_R_pt = Interp_Y_FD(0,0)
                if log_ind:
                    ax.fill_between([log10(No_R_pt - PARAMS.Yield_threshold),10],-10,10, facecolor='grey')
            x_points, y_points = np.zeros((Y2.shape[0]*Y2.shape[1]+1,n_phi,n_phi)), np.zeros((Y2.shape[0]*Y2.shape[1]+1,n_phi,n_phi))
            for ii in range(1,n_phi):
                for jj in range(1,n_phi):
                    if Interp_Y_FD(phi_1_values[ii],phi_2_values[jj])>PARAMS.Yield_threshold and phi_1_values[ii]/(phi_1_values[ii]+phi_2_values[jj])>0.5: # only care if we start from acceptable phi
                        x_points[0,ii,jj] = Interp_Y_FD(phi_1_values[ii],phi_2_values[jj])
                        y_points[0,ii,jj] = phi_1_values[ii]/(phi_1_values[ii]+phi_2_values[jj])
                        k = 1
                        for i in range(Y2.shape[0]): # doses
                            for j in range(Y2.shape[1]): # doses
                                x_points[k,ii,jj] = Y2[i,j,ii,jj]
                                y_points[k,ii,jj] = Res_A1[i,j,1,ii,jj]/(Res_A1[i,j,1,ii,jj]+Res_A2[i,j,1,ii,jj])
                                k = k+1
                        if not log_ind:
                            ax.scatter(x_points[0,ii,jj],y_points[0,ii,jj],color='k',marker='s')
                            ax.scatter(x_points[:,ii,jj],y_points[:,ii,jj],color= cmap(col/(n_phi^2)),marker='+')
                        if log_ind:
                            log_x_pts = [log10(-x_points[j,ii,jj]+No_R_pt) for j in range(x_points.shape[0])]
                            ax.scatter(log10(-x_points[0,ii,jj]+No_R_pt),y_points[0,ii,jj],color='k',marker='s')
                            ax.scatter(log_x_pts,y_points[:,ii,jj],color= cmap(col/(n_phi^2)),marker='+')
                        col = col + 1
        else:
            # Rv1, Rv2, Res_A1, Res_A2, dist1, dist2, Yield, Y2, SR1, SR2, phi_1_values, phi_2_values, Yield_vec, Diff, Interpolate_Y_FD = phi_plane_jumps(C_cont = Y_cont,n_phi=n_phi,n_d=n_doses,n_steps=n_steps)
            pp_jumps_2_out = phi_plane_jumps(C_cont = Y_cont,n_phi=n_phi,n_d=n_doses,n_steps=n_steps)
            Y2           = pp_jumps_2_out['Y2']
            Res_A1       = pp_jumps_2_out['Res_array_1']
            Res_A2       = pp_jumps_2_out['Res_array_2']
            phi_1_values = pp_jumps_2_out['phi_1_values']
            phi_2_values = pp_jumps_2_out['phi_2_values']
            if Interpolate_Y_FD is None:
                Interp_Y_FD = pp_jumps_2_out['Interpolate_Y_FD']
                No_R_pt = Interp_Y_FD(0,0)
                if log_ind:
                    ax.fill_between([log10(No_R_pt - PARAMS.Yield_threshold),10],-10,10, facecolor='grey')
            x_points, y_points = np.zeros((Y2.shape[0]*Y2.shape[1]+1,n_steps)), np.zeros((Y2.shape[0]*Y2.shape[1]+1,n_steps))
            for ii in phi_vector:
                x_points[0,ii] = Interp_Y_FD(phi_1_values[ii],phi_2_values[ii])
                y_points[0,ii] = phi_1_values[ii]/(phi_1_values[ii]+phi_2_values[ii])
                k = 1
                for i in range(Y2.shape[0]): # doses
                    for j in range(Y2.shape[1]): # doses
                        x_points[k,ii] = Y2[i,j,ii]
                        y_points[k,ii] = Res_A1[i,j,1,ii]/(Res_A1[i,j,1,ii]+Res_A2[i,j,1,ii])
                        if dose_descriptors:
                            Y2j = Y2[i,:,ii]
                            if len(np.argwhere(Y2j>50))>0:
                                j_min = min(np.argwhere(Y2j>50))
                                if j == j_min and j>0:
                                    if not log_ind:
                                        ax.scatter(x_points[k,ii],y_points[k,ii],color='g',marker='D') # contour dose
                                    if log_ind:
                                        ax.scatter(log10(-x_points[k,ii]+No_R_pt),y_points[k,ii],color='g',marker='D') # contour dose
                            if not log_ind:
                                if i == 0 or j == 0:
                                    ax.scatter(x_points[k,ii],y_points[k,ii],color='r',marker='D') # 0 dose
                                if i == Y2.shape[0]-1 or j == Y2.shape[1]-1:
                                    ax.scatter(x_points[k,ii],y_points[k,ii],color='b',marker='D') # full dose
                                if i == Y2.shape[0]-1 and j == Y2.shape[1]-1:
                                    ax.scatter(x_points[k,ii],y_points[k,ii],color='k',marker='D') # full dose both
                            if log_ind:
                                if i == 0 or j == 0:
                                    ax.scatter(log10(-x_points[k,ii]+No_R_pt),y_points[k,ii],color='r',marker='D') # 0 dose
                                if i == Y2.shape[0]-1 or j == Y2.shape[1]-1:
                                    ax.scatter(log10(-x_points[k,ii]+No_R_pt),y_points[k,ii],color='b',marker='D') # full dose
                                if i == Y2.shape[0]-1 and j == Y2.shape[1]-1:
                                    ax.scatter(log10(-x_points[k,ii]+No_R_pt),y_points[k,ii],color='k',marker='D') # full dose both
                        k = k+1
                ##
                if not log_ind:
                    ax.scatter(x_points[0,ii],y_points[0,ii],color='k',marker='s')
                    ax.scatter(x_points[:,ii],y_points[:,ii],color= cmap(col/(n_steps+1)),marker='+')
                if log_ind:
                    log_x_pts = [log10(-x_points[j,ii]+No_R_pt) for j in range(x_points.shape[0])]
                    ax.scatter(log10(-x_points[0,ii]+No_R_pt),y_points[0,ii],color='k',marker='s')
                    ax.scatter(log_x_pts,y_points[:,ii],color= cmap(col/(n_steps+1)),marker='+')
                col = col + 1
    else:
        marker_index = ['$0$','$1$','$2$','$3$','$4$','$5$','$6$','$7$','$8$','$9$','$10$','$11$','$12$','$13$','$14$','$15$','$16$','$17$','$18$','$19$','$20$']
        marker_index = marker_index[0:len(hlines)]
        ax.scatter(x_points[0],y_points[0],color='k',marker='s',s=50) # colorvec
        for k in range(1,len(x_points)):
            if colour_vec[k] > 0: # so that first year just has crosses
                ax.scatter(x_points[k],y_points[k],color=cmap(colour_vec[k]/(n_s-1)),marker=marker_index[int(bin_from[k])],s=50)
            else:
                ax.scatter(x_points[k],y_points[k],color=cmap(colour_vec[k]/(n_s-1)),marker='x',s=50)
        # find the route taken
        n_years = floor(len(x_points)/(len(hlines)-1))
        for k in range(1,n_years+1):
            k2 = (k-1)*(len(hlines)-1)
            x_pts = x_points[k2+1:k2+len(hlines)]
            y_pts = y_points[k2+1:k2+len(hlines)]
            x_pts = x_pts[x_pts>0]
            y_pts = y_pts[y_pts>0]
            ax.plot(x_pts,y_pts,color=cmap(colour_vec[k2]/(n_s-1)),label='Year % r' % (k-1))
        k_final = len(x_points) - len(hlines) + 1
        x_pts_f = max(x_points[k_final:])
        if x_pts_f>90:
            indices = np.zeros(n_years+1)
            ind = np.argwhere(x_points==x_pts_f)
            indices[0] = ind

            for kk in range(1,n_years):
                ind2 = int(bin_from[int(indices[kk-1])])
                k2 = (n_years - 1 - kk)*(len(hlines)-1) + 1
                # print(k2,int(indices[kk-1]),ind2)
                indices[kk] = k2 + ind2
            ##
            x_path = [x_points[int(jj)] for jj in indices]
            y_path = [y_points[int(jj)] for jj in indices]
            ax.plot(x_path,y_path,color='g',label='Path')

        ax.legend()
    ##


    # lower_lim = lower_lim - (100-upper_lim)-0.1 # min(lower_lim - (100-upper_lim)-0.1,PARAMS.Yield_threshold-1)
    # ax.set_xlim((lower_lim,100.1))
    upper_lim = np.amax(x_points)
    if x_lim is None:
        x_points[x_points<80] = 110
        lower_lim = np.amin(x_points)
        extra = (upper_lim-lower_lim)*0.15
        x_lo = lower_lim-extra
        x_hi = upper_lim+extra
    else:
        x_lo = x_lim[0]
        x_hi = x_lim[1]
    
    if not log_ind:
        ax.set_xlim((x_lo,x_hi))
    if log_ind:
        ax.set_xlim((log10(-upper_lim + No_R_pt)-0.5 , 1.5))

    ax.set_ylim((-0.05,1.05))
    
    if hlines is not None:
        for i in range(len(hlines)):
            ax.axhline(hlines[i],linestyle='--',color='k',alpha=0.5)
            if i < len(hlines)-1:
                ax.scatter((x_hi+No_R_pt)*0.5,0.5*(hlines[i+1]+hlines[i]),marker=marker_index[i],color='b')
    #
    ax.set_ylabel(r'Proportion of resistance to $F_1$',fontsize=15) #$\frac{\phi_1}{\phi_1+\phi_2}$
    if not log_ind:
        ax.set_xlabel('Yield Contour (FD)',fontsize=15)
        ax.invert_xaxis()
    if log_ind:
        ax.set_xlabel(r'Resistance Metric',fontsize=15) #$log_{10} (YC_{R0} - YC)$ (FD)
    # ax.grid()
    if title:
        if Y_cont is not None:
            ax.set_title('Yield contour, RF ratio; YC = %s' % round(Y_cont,3),fontsize=15)
        else:
            ax.set_title('Yield contour, RF ratio',fontsize=15)
    return None


# #----------------------------------------------------------------------------------------------      
def Y_t_space_contour_finder(Interpolate_Y_at_FD=None,res1=10**(-4),res2=10**(-4),n_d=4,t_vector=[0,0.1,0.5,0.9,1],n_seasons=2):
    n_seas = 1
    if Interpolate_Y_at_FD is None:
        Interpolate_Y_at_FD = Int_Y_FD(n_phi = 5)[0] # , Int_Y_at_FD
    ## Year i
    # generate up to n = 5*N points around the dose-space boundary ... linspaced on straight boundaries, and for n different t values on contour. Use interpolated fn to establish which on the straight boundary can be ignored
    # use these points to generate n new RF pairs
    n_bin = 1
    dose_vec = np.zeros((2,n_d))
    dose_vec[0,:] = [0.5*(i/(n_d-1)) for i in range(n_d)]
    dose_vec[1,:] = [0.5*(i/(n_d-1)) for i in range(n_d)]
    ##
    dose_used = -1*np.ones((2,len(t_vector)-1,n_seasons))
    res_used  =  1*np.ones((2,len(t_vector)-1,n_seasons))
    t_used, Y2_used, bin_from    = [-1*np.ones((len(t_vector)-1,n_seasons)) for i in range(3)]
    ##
    for seas_no in range(n_seasons):
        if seas_no>0:
            n_bin = len(t_vector)-1
        ##
        print('Season number: ',seas_no)
        ##
        contour   = -1*np.ones((2,n_d)) # needs mechanism for finding properly
        dose_gen, res_gen  = [-1*np.ones((2,n_bin,5*n_d-4)) for i in range(2)]
        ##
        t_gen, Y1_gen, Y2_gen = [-1*np.ones((n_bin,5*n_d-4)) for i in range(3)]
        ##
        for bin in range(n_bin):
            if seas_no==0:
                r1, r2 = res1, res2
            else:
                res_freq = res_used[:,:,seas_no-1]
                r1, r2 = res_freq[0,bin], res_freq[1,bin]
            ###
            # LTY, TY, FY, Yield_array, Res_array_1, Res_array_2, PRR_array, PRS_array, PSR_array, PSS_array, S1, S2, Innoc_array = master_loop_grid_of_tactics(n_d,n_seas,r1,r2,yield_stopper=PARAMS.Yield_threshold-0.5)
            grid_output = master_loop_grid_of_tactics(n_d,n_seas,r1,r2,yield_stopper=PARAMS.Yield_threshold-0.5)
            Yield_array = grid_output['Yield']
            Res_array_1 = grid_output['Res_array_1']
            Res_array_2 = grid_output['Res_array_2']

            k = 0
            for j in [0,-1]: # some sort of criteria to where poor doses not scanned
                for i in range(n_d):
                    if Yield_array[i,j,0]>PARAMS.Yield_threshold:
                        dose_gen[0,bin,k] = dose_vec[0,i]
                        dose_gen[1,bin,k] = dose_vec[1,j]
                        res_gen[0,bin,k]  = Res_array_1[i,j,1]
                        res_gen[1,bin,k]  = Res_array_2[i,j,1]
                        t_gen[bin,k]      = res_gen[0,bin,k]/(res_gen[0,bin,k]+res_gen[1,bin,k])
                        Y1_gen[bin,k]     = Yield_array[i,j,0]
                        Y2_gen[bin,k]     = Interpolate_Y_at_FD(Res_array_1[i,j,1],Res_array_2[i,j,1])
                        k = k+1
            for j in [0,-1]: # some sort of criteria to where poor doses not scanned?
                for i in range(1,n_d-1):
                    if Yield_array[j,i,0]>PARAMS.Yield_threshold:
                        dose_gen[0,bin,k] = dose_vec[0,j]
                        dose_gen[1,bin,k] = dose_vec[1,i]
                        res_gen[0,bin,k] = Res_array_1[j,i,1]
                        res_gen[1,bin,k] = Res_array_2[j,i,1]
                        t_gen[bin,k]     = res_gen[0,bin,k]/(res_gen[0,bin,k]+res_gen[1,bin,k])
                        Y1_gen[bin,k]    = Yield_array[j,i,0]
                        Y2_gen[bin,k]    = Interpolate_Y_at_FD(Res_array_1[j,i,1],Res_array_2[j,i,1])
                        k = k+1
            # define interpolated function
            Interpolate_Y = interpolate(Yield_array[:,:,0],Yield_array.shape[0],n_d,'linear')[0] # , Int_Y
            def equations(y):
                return Interpolate_Y(x_val,y)-PARAMS.Yield_threshold
            ##
            ig = 0.2
            ##
            for i in range(n_d):
                x_val = i/(n_d-1) # in [0,1]
                y_val = fsolve(equations,ig) # in [0,1]          # given x, find y such that are on 95% contour
                if y_val == ig:
                    y_val = fsolve(equations,0.5)
                if y_val<0 or y_val>1 or y_val==0.5:
                    warnings.warn('Warning! y_val error: y_val = ',y_val,i) # print # was None with this message commented out
                else:
                    contour[0,i] = 0.5*x_val
                    contour[1,i] = 0.5*y_val
                    # print('Contour used: ',contour[:,i])
                ##
                    one_t_output =  master_loop_one_tactic([contour[0,i]],[contour[1,i]],[contour[0,i]],[contour[1,i]],r1,r2,yield_stopper=PARAMS.Yield_threshold-0.5)
                    Yield_vec = one_t_output['Yield_vec']
                    Res_vec_1 = one_t_output['Res_vec_1']
                    Res_vec_2 = one_t_output['Res_vec_2']
                    ##
                    dose_gen[:,bin,k] = contour[:,i]
                    res_gen[0,bin,k]  = Res_vec_1[1]
                    res_gen[1,bin,k]  = Res_vec_2[1]
                    t_gen[bin,k]      = res_gen[0,bin,k]/(res_gen[0,bin,k]+res_gen[1,bin,k])
                    Y1_gen[bin,k]     = Yield_vec[0]
                    Y2_gen[bin,k]     = Interpolate_Y_at_FD(Res_vec_1[1],Res_vec_2[1])
                    k = k+1
        ###
        # print(contour)
        # print(Y1_gen,'y1gen')
        # print(Y2_gen,'y2gen')
        # once each year
        # create k (or k*n) different bins in which we select the lowest Y2 value for each - probably more bins near t = 1
        for kk in range(len(t_vector)-1): # for t in the relevant bin, find minimum Y2 contour
            index = np.argwhere((t_gen<t_vector[kk+1])&(t_gen>t_vector[kk]))
            ##
            Y_top = 0
            if index.shape[0]>0:
                Y2_top = [Y2_gen[index[i][0],index[i][1]] for i in range(index.shape[0])]
                for i in range(index.shape[0]):
                    if Y1_gen[index[i][0],index[i][1]]>PARAMS.Yield_threshold: # acceptable dose
                        if Y2_gen[index[i][0],index[i][1]]>Y_top:
                            Y_top = Y2_gen[index[i][0],index[i][1]]
                i2 = np.argwhere(Y2_top==Y_top)[0][0]
                ##
                dose_used[:,kk,seas_no] = dose_gen[:,index[i2][0],index[i2][1]]
                res_used[:,kk,seas_no]  = res_gen[:,index[i2][0],index[i2][1]]
                t_used[kk,seas_no]      = t_gen[index[i2][0],index[i2][1]]
                Y2_used[kk,seas_no]     = Y2_gen[index[i2][0],index[i2][1]]
                bin_from[kk,seas_no]    = index[i2][0]
        # save the relevant dose values and RF values
    # print(Y2_used)
    # print(t_used)
    # print(res_used)
    # print(dose_used)
    # print(bin_from)
    return t_used, Y2_used, res_used, dose_used, bin_from




# #----------------------------------------------------------------------------------------------      
def FD_space(Yield_out,bottom_phi):
    n_p = Yield_out.shape[1]
    phi_vec = np.linspace(bottom_phi,0,n_p)
    phi_vec_rr = np.linspace(2*bottom_phi,0,2*n_p)
    Y2_array = 50*np.ones((2*n_p,n_p,n_p))
    F_array   = np.zeros((2*n_p,n_p,n_p))
    for i in range(2*n_p):
        for j in range(n_p):
            for k in range(n_p):
                prr = 10**(phi_vec_rr[i])
                prs = 10**(phi_vec[j])
                psr = 10**(phi_vec[k])
                pss = 1 - prr - prs - psr
                if pss>0:
                    Y2_array[i,j,k] = Yield_out[i,j,k,-1,-1]
                if Y2_array[i,j,k]>PARAMS.Yield_threshold:
                    F_array[i,j,k] = 1
    x   = phi_vec_rr
    y,z = [phi_vec for i in range(2)]
    Int_Y2   = RegularGridInterpolator((x,y,z),Y2_array)
    dictionary = {'F_array': F_array, 'Y2_array': Y2_array, 'Int_Y2': Int_Y2}
    return dictionary






##----------------------------------------------------------
def C_contours(n_p=18,n_c=18,p_bottom=10**(-4),p_top=1-10**(-4),c_bottom=10**(-8),c_top=0.5,dose=0.5,C_vec=None,p=None,phi_vec=None,Y_bottom=PARAMS.Yield_threshold-0.1,surface_plot=None):
    ##----------------------------------------------------------
    if C_vec is None:
        C_vec = np.linspace(c_bottom,c_top,n_c)
    if p is None:
        p = np.linspace(p_bottom,p_top,n_p)
    ##----------------------------------------------------------
    Yield_vec = Y_bottom*np.ones((n_p,len(C_vec)))
    phi_here = np.zeros((n_p,len(C_vec)))
    # phi_here1 = np.zeros((n_p,len(C_vec)))
    # phi_here2 = np.zeros((n_p,len(C_vec)))
    bound_L = np.zeros((n_p,len(C_vec)))
    bound_R = np.zeros((n_p,len(C_vec)))
    ##----------------------------------------------------------
    for i in range(n_p):
        for j in range(len(C_vec)):
            phi_here[i,j] = (C_vec[j]/(p[i]*(1-p[i])))**(0.5)
            bound_L[i,j] = - phi_here[i,j] + 1/(1-p[i])
            bound_R[i,j] = - phi_here[i,j] + 1/(p[i])
            if phi_here[i,j]<min(1/p[i],1/(1-p[i])):
                output = master_loop_one_tactic([dose],[dose],[dose],[dose],p[i]*phi_here[i,j],(1-p[i])*phi_here[i,j]) # but full dose might be 0.5
                Yield_vec[i,j] = output['Yield_vec']

    if phi_vec is not None:
        p_scat = [((phi_vec[0,i]/(phi_vec[0,i]+phi_vec[1,i])) - p_bottom)/(p_top-p_bottom) for i in range(len(phi_vec[0,:]))]
        C_scat = [(phi_vec[0,i]*phi_vec[1,i] - c_bottom)/(c_top-c_bottom) for i in range(len(phi_vec[0,:]))] # scaled to match axis scale - won't be perfect but won't be far off
        Con_at = {'Con_mats':(Yield_vec,bound_L,bound_R),'Con_inline':('inline',None,None),'Con_levels':([85,90,91,92,93,94,95,96,97,98,99],[-10**(30),0],[-10**(30),0]),'Contourf': (None,'k','k')}
        Overlay_plotter(Con_attr=Con_at,Scat_attr = {'x_scat': p_scat,'y_scat':C_scat}, figure_attributes={'x_lab':'p','y_lab':'C','xtick': p,'ytick': C_vec,'title':'Yield, p, C'})
    else:
        Overlay_plotter(Con_attr={'Con_mats':(Yield_vec,bound_L,bound_R),'Con_inline':('inline',None,None),'Con_levels':([85,90,91,92,93,94,95,96,97,98,99],[-10**(30),0],[-10**(30),0]),'Contourf': (None,'k','k')},figure_attributes={'x_lab':'p','y_lab':'C','xtick': p,'ytick': C_vec,'title':'Yield, p, C'})

    # for i in range(n_p):
    #     phi_here1[i,:] = phi_here[i,:]*(p[i]) # phi1
    #     phi_here2[i,:] = phi_here[i,:]*(1-p[i]) # phi2
    # Overlay_plotter(Col_mat= Yield_vec,Col_label='Yield',Col_bds_vec=(PARAMS.Yield_threshold,),title='Yield, p, C')
    # Overlay_plotter(Con_mats=(phi_here),Con_levels=([0,0.2,0.4,0.6,0.8,1]),Con_inline=('inline'),x_lab='p',y_lab='C',xtick= p,ytick= C_vec,title='phi_1 + phi_2, p, C')
    # Overlay_plotter(Con_mats=(phi_here1,phi_here2),Con_levels=([0,0.1,0.2,0.3,0.4,0.5],[0,0.1,0.2,0.3,0.4,0.5]),Con_inline=('inline','inline'),x_lab='p',y_lab='C',xtick= p,ytick= C_vec,title='phi_1, phi_2, p, C')

    if surface_plot is not None:    
        X, Y = np.meshgrid(np.linspace(0,1,n_c),np.linspace(0,1,n_p))
        fig5 = plt.figure()
        ax5 = fig5.gca(projection='3d')
        ax5.plot_wireframe(X,Y,Yield_vec)
    return None