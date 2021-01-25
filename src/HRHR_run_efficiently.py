import sys
sys.path.insert(0, '/Functions_and_plotting/')
##
import matplotlib.pyplot as plt
import numpy as np
from math import ceil, log, log10, floor
from scipy.interpolate import interp1d
from matplotlib.ticker import MaxNLocator
from scipy.optimize import fsolve
##
from Functions_and_plotting.plotter_HRHR    import Overlay_plotter, plot_disease_dynamics, bar_plotter, line_or_scatter_plotter # Equal_dose,
from Functions_and_plotting.functions_HRHR  import interpolate, solv_ode, master_loop_one_tactic, master_loop_grid_of_tactics, Z_metric, Dose_tuplet_extractor, fngcide
from Functions_and_plotting.parameters_HRHR import params
#----------------------------------------------------------------------------------------------
Y_levels = [95,96,97,98,99]
Y_levels_sparse = [95,96,97]
Z_levels = [0.3,0.4,0.5,0.6,0.7]
Z_levels_sparse = [0.5]
LTY_levels = [10,12,15]
SR_levels  = [1.2,1.6,2,2.8]
SR_levels2  = [1.6,1.9,2.2,2.5,2.8]
# #----------------------------------------------------------------------------------------------
### for plot_Y_SR_RF_by_year
i_vector = [0.25,0.5,0.75,1] #[0.4,0.5,0.6,0.8,0.9,1]
j_vector = [0.25,0.5,0.75,1] #[0.4,0.5,0.6,0.8,0.9,1]
### for dynamics_by_SR_RF
r_vec = [10**(-6),10**(-3),0.01,0.02,0.04,0.07] #[10**(-6),10**(-3),10**(-2),0.02,0.03,0.04,0.05,0.06,0.07,0.1,0.4] # needs to be length n_here
n_here= len(r_vec)
n_d = 20 # does it have to be n_here??? don't think so
d_vec = np.linspace(0,2,n_d)
# d_vec = [0,0.05,0.1,0.15,0.3,0.7,0.85,0.9,0.95,1] # other way to specify d_vec, but then the plotter is less effective
# see d1, r1... below
# #----------------------------------------------------------------------------------------------
y_thresh = 70
r1 = 10**(-5)
r2 = 10**(-5)
n_doses = 5
n_seasons = 35
strat = 'mix'
output_size = 101 # for interpolation
jet_cmap = plt.get_cmap('jet')

# should update so that the exponential decays are just functions
#  write as a class

# #----------------------------------------------------------------------------------------------
# plots
delt_quad_for_FYR     = False
DR_for_pres           = False
plot_Y_SR_RF_by_year  = False   # x: year; y: Y/SR/RF; for combination of doses specified by i_vector and j_vector. Can also do log(RF) on y-axis
interpolate_and_plots = False    # x/y: doses; z: LTY/TY. Plots colour/contour/overlay/wireframe plots, some using interpolated function
Y_SR_RF_by_equal_dose = False   # x: dose; y: Y/SR/RF. For multiple years, d1=d2=dose
bar_plots_RFs         = False    # fig 6 poster; x: time; y: disease dynamics, for doses, RFs specified by r_vec, d_vec.
plot_dis_dyn          = False    # disease dynamics
plot_SR_equal_dose    = False    # fig 5 poster Also colour/wireframe/contour plots with dose and RF on x/y, SR and Y on z. Also Y/SR with x: dose, varying RFs.
plot_SR_Z_by_year     = False    # x/y: doses, z: SR/Z. Contour plots of SR/Z for different years
contour_plots         = False    # problem   #1 # x/y: doses, z: various
colour_plots          = False    # x/y: doses, z: various
LTY_plot              = True     # x/y: doses, z: LTY
FY_plot               = True     # x/y: doses, z: FY
breakdown_year        = True
# #----------------------------------------------------------------------------------------------
# #----------------------------------------------------------------------------------------------
# ignore these since they will be made True if needed
overlay_plots         = False     # default
run_simulation        = False     # run simulation which loops over grid of doses
dynamics_by_SR_RF     = False     #1 #1 # Runs own simulations
##
chapt_3_FYR = False
if chapt_3_FYR:
        plot_SR_equal_dose = True
        delt_quad_for_FYR = True
        FY_plot = True
if LTY_plot or FY_plot or breakdown_year:
        overlay_plots = True    
if plot_Y_SR_RF_by_year or interpolate_and_plots or Y_SR_RF_by_equal_dose or plot_SR_Z_by_year or contour_plots or colour_plots or overlay_plots:
        run_simulation = True
if bar_plots_RFs or plot_SR_equal_dose:
        dynamics_by_SR_RF = True
# #----------------------------------------------------------------------------------------------
# delta quadratic for FYR
if delt_quad_for_FYR:
        x = np.linspace(0,1,100)
        y = np.zeros(100)

        for i in range(len(x)):
                y[i] = x[i]*(1-x[i])

        fig = plt.figure(figsize=(9,6))
        ax = fig.add_subplot(111)
        ax.plot(x,y)
        ax.set_ylabel(r'$\frac{1}{\rho}$ln($SR$)',fontsize=16)
        ax.set_xlabel(r'$\delta$;'+' reduction in growth rate',fontsize=16)
        ax.text(0.87,-0.12,'Low dose',horizontalalignment='center',verticalalignment='bottom',fontsize=14,transform=ax.transAxes)
        ax.text(0.13,-0.12,'High dose',horizontalalignment='center',verticalalignment='bottom',fontsize=14,transform=ax.transAxes)
        ax.grid()
# #----------------------------------------------------------------------------------------------
# Dose Response for Presentation
if DR_for_pres:
        x = np.linspace(0,1,100)
        y = np.zeros(100)

        for i in range(len(x)):
                y[i] = 100*fngcide(params.omega_1,params.theta_1,x[i])
                # y[i] = 100*fngcide(0.48,9.9,x[i])


        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(x,y)
        ax.set_ylabel('Effective Infection Rate (%)',fontsize=16)
        ax.set_xlabel('Fraction of Maximum Legal Dose',fontsize=16)
        # ax.set_title('Dose Response: Chlorothalonil',fontsize=16)
        ax.set_title('Dose Response: Pyraclostrobin',fontsize=16)
        ax.grid()
        ax.set_xlim((0,1))
        ax.set_ylim((0,100))
# #----------------------------------------------------------------------------------------------
##
if run_simulation:
        output = master_loop_grid_of_tactics(n_doses,n_seasons,r1,r2,y_thresh)
        LTY = output['LTY']
        TY = output['TY']
        FY = output['FY']
        Yield = output['Yield']
        Res_array_1 = output['Res_array_1']
        Res_array_2 = output['Res_array_2']
        prr = output['PRR_array']
        prs = output['PRS_array']
        psr = output['PSR_array']
        pss = output['PSS_array']
        Selection_array_1 = output['Selection_array_1']
        Selection_array_2 = output['Selection_array_2']
        Innoc_array = output['Innoc_array']
        t1 = output['t_vec']
        Z = Z_metric(Selection_array_1[:,:,1],Selection_array_2[:,:,1])
        FYY = Yield[:,:,0]
# #----------------------------------------------------------------------------------------------
if plot_Y_SR_RF_by_year:
        logRes1 = np.zeros((Res_array_1.shape[0],Res_array_1.shape[1],Res_array_1.shape[2]))
        for i in range(Res_array_1.shape[0]):
                for j in range(Res_array_1.shape[1]):
                        for k in range(Res_array_1.shape[2]):
                                logRes1[i,j,k] = log(Res_array_1[i,j,k])

        SR_tup1,SR_tup2,R_tup1,R_tup2,Y_tup,L_tup,C_tup = Dose_tuplet_extractor(Selection_array_1,Selection_array_2,Res_array_1,Res_array_2,Yield,i_vector,j_vector,n_doses)
        # SR_tup1,SR_tup2,LR_tup1,R_tup2,Y_tup,L_tup,C_tup = Dose_tuplet_extractor(Selection_array_1,Selection_array_2,logRes1,Res_array_2,Yield,i_vector,j_vector,n_doses)

        subplot_1 = {'line_y_data': [Y_tup[i] for i in range(len(Y_tup))],
        'xlab': 'Season Number',
        'ylab': r'Yield ($\%$ of DF)',
        'xlim': (0,len(Y_tup[0])),
        'ylim': (np.amin(Y_tup)-1,100),
        'grey_x_box_pos': [-1,2*len(Y_tup[0])],
        'grey_y_box_start': -1,
        'grey_y_box_end': params.Yield_threshold,
        'line_label': L_tup,
        'line_colour': C_tup}

        subplot_2 = {'line_x_data': [np.linspace(0,len(R_tup1[0])-1,len(R_tup1[0]))]*len(R_tup1),
        'line_y_data': [R_tup1[i] for i in range(len(R_tup1))],
        'xlim': (0,len(R_tup1[0])-1),
        'xlab': 'Season Number',
        'ylab': r'RF ($F_1$)',
        'line_colour': C_tup}

        figure_attributes = {'legend': 1, 'leg_pos': 'upper right', 'horiz': 1, 'fs': (18,5), 'label_fontsize': 16}

        line_or_scatter_plotter(plots = [subplot_1,subplot_2],figure_attributes = figure_attributes)
# #----------------------------------------------------------------------------------------------
####
if interpolate_and_plots:

        Interpolate, Int = interpolate(LTY,n_doses,output_size,'linear')

        Col_attr_1 = {'Col_mat': LTY, 'Col_label': 'LY', 'Col_bds_vec': (0,)}
        Col_attr_2 = {'Col_mat': Int, 'Col_label': 'LY', 'Col_bds_vec': (0,)}
        Col_attr_3 = {'Con_mats': (Int), 'Con_inline': ('inline')}
        Col_attr_4 = {'Col_mat': Int, 'Col_label': 'Yield', 'Col_bds_vec': (0,)} 
        Con_attr_4 = {'Con_mats': (Z,FYY), 'Con_levels': (Z_levels_sparse,Y_levels), 'Con_inline': (None,'inline')}
        Overlay_plotter(Col_attr=Col_attr_1, figure_attributes={'title': 'not_interp'})
        Overlay_plotter(Col_attr=Col_attr_2, figure_attributes = {'title': '!! - warning interp LY'})
        Overlay_plotter(Col_attr=Col_attr_3, figure_attributes = {'title': '!! - warning interp LY'})
        Overlay_plotter(Col_attr=Col_attr_4,Con_attr=Con_attr_4, figure_attributes = {'title': 'LTY, Mixture'})
        ####
        X, Y = np.meshgrid(np.linspace(0,1,output_size),np.linspace(0,1,output_size))
        fig2 = plt.figure()
        ax2 = fig2.gca(projection='3d')
        surf2 = ax2.plot_wireframe(X,Y,Int)
        ax2.set_title('LTY')
        ####
        Interpolate, Int_TY = interpolate(TY,n_doses,output_size,'linear') # cubic gives neg values

        # Overlay_plotter(TY,'TY',(17,),figure_attributes={'title':'not_interp'})
        # # Colour_plotter(Int_TY,'TY',(17,),title='!! - warning interp')
        # Overlay_plotter(Con_mats= (Int_TY), figure_attributes = {'title': '!! - warning interp TY'}, Con_inline= ('inline'))

        fig2a = plt.figure()
        ax2a = fig2a.gca(projection='3d')
        surf2a = ax2a.plot_wireframe(X,Y,Int_TY)
        ax2a.set_title('TY')
# #----------------------------------------------------------------------------------------------
### Yield/SR/RF for d1=d2 varying and years increasing
if Y_SR_RF_by_equal_dose:
        # season_step = 2
        # Equal_dose(Yield,Selection_array_1,logRes1,season_step,n_doses,n_seasons)
        # Equal_dose(Yield,Selection_array_1,Res_array_1,season_step,n_doses,n_seasons)

        SR1 = np.zeros((n_doses,n_seasons))
        Y   = np.zeros((n_doses,n_seasons))
        RF  = np.zeros((n_doses,n_seasons))
        for ii in range(n_doses):
                for jj in range(n_seasons):
                        SR1[ii,jj] = Selection_array_1[ii,ii,jj]
                        Y[ii,jj]   = Yield[ii,ii,jj]
                        RF[ii,jj]   = Res_array_1[ii,ii,jj]

        subplot_1 = {'line_y_data': [Y[:,i] for i in range(Y.shape[1])],
        'line_x_data': [np.linspace(0,1,Y.shape[0]) for i in range(Y.shape[1])],
        'ylab': 'Yield',
        'line_label': ['Year %s' % jj for jj in range(Y.shape[1])],
        'line_colour': [jet_cmap(jj/(Y.shape[1]-1)) for jj in range(Y.shape[1])]}

        subplot_2 = {'line_y_data': [SR1[:,i] for i in range(SR1.shape[1])],
        'line_x_data': [np.linspace(0,1,SR1.shape[0]) for i in range(SR1.shape[1])],
        'ylab': 'SR',
        'line_colour': [jet_cmap(jj/(SR1.shape[1]-1)) for jj in range(SR1.shape[1])]}

        subplot_3 = {'line_y_data': [RF[:,i] for i in range(RF.shape[1])],
        'line_x_data': [np.linspace(0,1,RF.shape[0]) for i in range(RF.shape[1])],
        'ylab': 'RF',
        'line_colour': [jet_cmap(jj/(RF.shape[1]-1)) for jj in range(RF.shape[1])]}

        figure_attributes = {'legend': 1, 'leg_pos': 'upper right'}

        line_or_scatter_plotter(plots = [subplot_1,subplot_2,subplot_3],figure_attributes = figure_attributes)
# #----------------------------------------------------------------------------------------------
if dynamics_by_SR_RF:
        SR_vec = np.zeros((len(d_vec),len(r_vec)))
        Y_vec = np.zeros((len(d_vec),len(r_vec)))
        k = 0
        no_years = 2
        Sol=np.zeros((ceil((params.T_GS87-params.T_emerge)/params.dt),params.no_variables,no_years,len(d_vec),len(r_vec)))
        t = np.zeros(ceil((params.T_GS87-params.T_emerge)/params.dt))
        for i in range(len(d_vec)):
                for j in range(len(r_vec)):
                        one_t_output = master_loop_one_tactic(d_vec[i]*np.ones(no_years),d_vec[i]*np.ones(no_years),d_vec[i]*np.ones(no_years),d_vec[i]*np.ones(no_years),r_vec[j],r_vec[j],y_thresh)
                        SR_vec[i,j]    = one_t_output['Selection_vec_1'][1]
                        Y_vec[i,j]     = one_t_output['Yield_vec'][0]
                        Sol[:,:,:,i,j] = one_t_output['Sol_array']
                        tt             = one_t_output['t_vec']

        #### 
if bar_plots_RFs:
        rr = [0,1,2,2,4]
        dd = [4,4,3,4,4]
        
        # Sol[t,variable,year,dose,res_freq]
# issues here
        scaled_sol = np.zeros((len(tt),params.no_variables,len(rr)))
        dydt = np.zeros((len(rr),params.no_variables,len(tt)))
        scaled_dydt = np.zeros((len(rr),params.no_variables,len(tt)))
        scaled_dydt2 = np.zeros((len(rr),params.no_variables,len(tt)))
        for jj in range(min(len(rr),len(dd))):
                for ii in range(len(tt)):
                        # print(Sol.shape)
                        # print(dydt.shape)
                        dydt[jj,:,ii] = solv_ode(tt[ii],Sol[ii,:,0,dd[jj],rr[jj]])
        
        # really should look to scale to find per capita growth rates

        for i in range(len(rr)):
                for k in [params.ER_ind,params.IR_ind]:
                        scaled_dydt[i,k,:] = (1/(r_vec[rr[i]]**2))*dydt[i,k,:]
                for j in [params.ERS_ind,params.ESR_ind,params.IRS_ind,params.ISR_ind]:
                        scaled_sol[:,j,i] = (1/(r_vec[rr[i]]))*Sol[:,j,0,dd[i],rr[i]]
                        scaled_dydt[i,j,:] = (1/(r_vec[rr[i]]))*dydt[i,j,:]
                        scaled_dydt2[i,j,:] = (1/(r_vec[rr[i]]))*dydt[i,j,:]
                for l in [params.ES_ind,params.IS_ind]:
                        scaled_sol[:,l,i] = Sol[:,l,0,dd[i],rr[i]]
                        scaled_dydt[i,l,:] = dydt[i,l,:]
                        scaled_dydt2[i,l,:] = dydt[i,l,:]
        
        # for i in range(len(rr)):
        # for i in range(len(rr)):
        #         plot_disease_dynamics(tt, scaled_sol[:,:,i],title = 'Scaled dyn., RF = %s, D = %s' % (r_vec[rr[i]],d_vec[dd[i]]),no_R='on')
        # for i in range(len(rr)):
        #         plot_disease_dynamics(t, np.transpose(dydt[i,:,:]),title = 'Derivatives, RF = %s, D = %s' % (r_vec[rr[i]],d_vec[dd[i]]))
        # for i in range(len(rr)):
        #         plot_disease_dynamics(t, np.transpose(scaled_dydt[i,:,:]),title = 'Scaled Deriv., RF = %s, D = %s' % (r_vec[rr[i]],d_vec[dd[i]]))
        # for i in range(len(rr)):
        #         plot_disease_dynamics(tt, np.transpose(scaled_dydt2[i,:,:]),title = 'Scaled Deriv., RF = %s, D = %s' % (r_vec[rr[i]],d_vec[dd[i]]),no_R='on')
        for i in range(len(rr)):
                bar_plotter(tt, Sol[:,:,0,dd[i],rr[i]],bar3=1,sub1=1) # title = 'RF = %s' % r_vec[rr[i]], # title = 'RF = %s, D = %s' % (r_vec[rr[i]],d_vec[dd[i]])

if plot_dis_dyn:
        d_vec1 = np.ones(1)*0.5
        d_vec2 = np.ones(1)*0.5

        one_output = master_loop_one_tactic(d_vec1,d_vec2,d_vec1,d_vec2,res_prop_1=10**(-1),res_prop_2=10**(-2))
        plot_disease_dynamics(one_output['t_vec'], one_output['Sol_array'][:,:,0])

if plot_SR_equal_dose:
        # SR_col_attr = {'Col_mat': SR_vec,'Col_label':'SR','Col_bds_vec': (1,)}
        # Overlay_plotter(Col_attr=SR_col_attr,figure_attributes = {'title': 'SR','x_lab': 'Dose','y_lab': 'RF'})
        
        Y_col_attr = {'Col_mat': Y_vec,'Col_label':'Y','Col_bds_vec': (90,)}
        Overlay_plotter(Col_attr= Y_col_attr,figure_attributes={'x_lab':'Dose','y_lab':'RF','ytick': r_vec,'title':'Y'})


        SR_con_attr = {'Con_mats': (SR_vec),'Con_inline':('inline')}
        Overlay_plotter(Con_attr=SR_con_attr,figure_attributes={'x_lab':'Dose','y_lab':'RF','ytick':r_vec,'title':'SR'})
        
        Y_con_attr = {'Con_mats': (Y_vec),'Con_inline':('inline')}
        Overlay_plotter(Con_attr=Y_con_attr,figure_attributes={'x_lab':'Dose','y_lab':'RF','ytick':r_vec,'title':'Y'})
        
        Con_attrs = {'Con_mats': (SR_vec,Y_vec),'Con_levels':(SR_levels2,[params.Yield_threshold]),'Con_inline':('inline',None)}
        Overlay_plotter(Con_attr=Con_attrs,figure_attributes={'x_lab':'Dose','y_lab':'RF','ytick':r_vec,'title':'Y'})
        
        # Con_attrs2 = {'Con_mats':(Y_vec),'Con_levels':([params.Yield_threshold]),'Con_inline':(None)}
        # Col_attrs2 = {'Col_mat': SR_vec,'Col_label':'SR','Col_bds_vec':(1,)}
        # Overlay_plotter(Con_attr=Con_attrs2,Col_attr=Col_attrs2,figure_attributes={'title':'SR','x_lab':'Dose','y_lab':'RF','ytick':r_vec})

        
        
        int_Y = interp1d(d_vec,Y_vec[:,0],kind='linear')
        def inter(d):
                res = int_Y(d) - params.Yield_threshold
                return res
        dose_95 = fsolve(inter,0.2)[0]

        lower_limit = dose_95

        subplot_1 = {'line_x_data': [d_vec]*n_here,
        'line_y_data': [SR_vec[:,jj] for jj in range(n_here)],
        'line_label': ['RF %s' % r_vec[jj] for jj in range(n_here)],
        'ylab': 'SR',
        'line_colour': [jet_cmap(jj/(n_here-1)) for jj in range(n_here)],
        'xlim': (0,max(d_vec)),
        'ylim': (1/1.05,1.05*np.amax(SR_vec)),
        'grey_x_box_pos': [-1,dose_95],
        'grey_y_box_start': -1,
        'grey_y_box_end': 100}

        subplot_2 = {'line_x_data': [d_vec]*n_here,
        'line_y_data': [Y_vec[:,jj] for jj in range(n_here)],
        'xlim': (0,max(d_vec)),
        'ylim': (np.amin(Y_vec)-1,100),
        'ylab': 'Yield',
        'line_colour': [jet_cmap(jj/(n_here-1)) for jj in range(n_here)],
        'grey_x_box_pos': [-1,2],
        'grey_y_box_start': -1,
        'grey_y_box_end': params.Yield_threshold}

        
        figure_attributes = {'legend': 1, 'leg_pos': 'lower right', 'fs': (10,7), 'share_x': True, 'label_fontsize': 16, 'horiz': 0, 'xlab': 'Dose'}

        line_or_scatter_plotter(plots = [subplot_1,subplot_2],figure_attributes = figure_attributes)

        X, Y = np.meshgrid(np.linspace(0,1,n_here),np.linspace(0,1,n_d))
        fig5 = plt.figure()
        ax5 = fig5.gca(projection='3d')
        surf5 = ax5.plot_wireframe(X,Y,SR_vec)
# #----------------------------------------------------------------------------------------------
#### contour plots of selection array and yield as year varies
if plot_SR_Z_by_year:
        for i in [1,6,11,16]:
                Z = Z_metric(Selection_array_1[:,:,i],Selection_array_2[:,:,i])
                Con_attrib = {'Con_mats': (Selection_array_1[:,:,i],Selection_array_2[:,:,i]),'Con_levels':(SR_levels,SR_levels),'Con_inline':('inline','inline')}
                Overlay_plotter(Con_attr=Con_attrib,figure_attributes={'title':'SR Year % s' % i})
                CA = {'Con_mats': (Z),'Con_levels': (Z_levels),'Con_inline':('inline')}
                Overlay_plotter(Con_attr=CA,figure_attributes={'title':'Z metric, Year % s' % i})
# # #----------------------------------------------------------------------------------------------
if contour_plots:
        CA = {'Con_mats':(LTY),'Con_inline':('Lifetime Yield')}
        Overlay_plotter(Con_attr=CA,figure_attributes={'title':'LTY'})

        CA = {'Con_mats':(LTY,FYY),'Con_levels':(LTY_levels,Y_levels),'Con_inline':('inline','FYY')}
        Overlay_plotter(Con_attr=CA,figure_attributes={'title':'LTY/Y'})

        CA = {'Con_mats':(Z,FYY),'Con_levels':(Z_levels_sparse,Y_levels_sparse),'Con_inline':(None,'inline')}
        Overlay_plotter(Con_attr=CA,figure_attributes={'title':'Y/Z'})

        CA = {'Con_mats':(FYY,Selection_array_1[:,:,1],Selection_array_2[:,:,1]),'Con_levels':(Y_levels_sparse,SR_levels,SR_levels),'Con_inline':('inline',None,'SR')}
        Overlay_plotter(Con_attr=CA,figure_attributes={'title':'Y/SR1/SR2'})

        CA = {'Con_mats':(FYY,Z,LTY),'Con_levels':(Y_levels,Z_levels_sparse,None),'Con_inline':('inline',None,'inline')}
        Overlay_plotter(Con_attr=CA,figure_attributes={'title':'LTY/Y/Z'})

# #----------------------------------------------------------------------------------------------
if colour_plots:
        CA = {'Col_mat':FYY,'Col_label':'Yield','Col_bds_vec': (90,)}
        Overlay_plotter(Col_attr=CA,figure_attributes={'title':'First Year Yield, Mixture'})
        CA = {'Col_mat':LTY,'Col_label':'Lifetime Yield','Col_bds_vec':(0,)}
        Overlay_plotter(Col_attr=CA,figure_attributes={'title':'LTY, Mixture'})
        CA = {'Col_mat':FY,'Col_label':'Failure Year','Col_bds_vec':(0,)}
        Overlay_plotter(Col_attr=CA,figure_attributes={'title':'FY, Mixture'})
        CA = {'Col_mat':TY,'Col_label':'Total Yield','Col_bds_vec':(0,)}
        Overlay_plotter(Col_attr=CA,figure_attributes={'title':'TY, Mixture'})
# # #----------------------------------------------------------------------------------------------
if overlay_plots:
        tick_vec = [0, 0.4,0.8,1.2,1.6,2] # or False
        alpha_use = 1
        if LTY_plot:
                print(LTY)
                LTY_min = np.amin(LTY[LTY>0])
                Cl_A = {'Col_mat': LTY,'Col_label': 'Lifetime Yield','Col_bds_vec': (0,LTY_min,),'alpha': alpha_use}
                Cn_A = {'Con_mats': (FYY,Z),'Con_inline':('inline',None),'Con_levels': ([95],[0.5])} # 'Con_levels': None,
                fig_A = {'xtick': tick_vec, 'ytick': tick_vec}
                Overlay_plotter(Col_attr=Cl_A, Con_attr=Cn_A,figure_attributes=fig_A)
                ##
                Cn_A = {'Con_mats': (Z),'Con_inline':(None),'Con_levels': ([0.5])} # 'Con_levels': None,
                Overlay_plotter(Col_attr=Cl_A, Con_attr=Cn_A,figure_attributes=fig_A)
        if FY_plot:
                FY_max = np.amax(FY)+1
                FY_min = np.amin(FY[FY>0])
                C2_A = {'Col_mat': FY,'Col_label': 'Failure Year','Col_bds_vec': (0,FY_min,),'alpha': alpha_use} # np.concatenate(0,np.linspace(FY_min,FY_max,FY_max-FY_min+1))}

                fig_A = {'xtick': tick_vec, 'ytick': tick_vec}
                Overlay_plotter(Col_attr=C2_A,figure_attributes=fig_A) # Con_attr=Cn_A) # Col_attr=Cl_A)
                ##
                Cn_A = {'Con_mats': (Z),'Con_inline':(None),'Con_levels': ([0.5])} # 'Con_levels': None,
                Overlay_plotter(Col_attr=C2_A, Con_attr=Cn_A,figure_attributes=fig_A)
                print(FY)

        
        if breakdown_year:
                
                dose_sum = range(1,2*n_doses-1)
                print(dose_sum)

                ratio = np.nan*np.zeros((n_doses,len(dose_sum)))
                fail  = np.nan*np.zeros((n_doses,len(dose_sum)))

                
                for i in range(len(dose_sum)):
                        for d1 in range(n_doses):
                                d1 = d1
                                d2 = dose_sum[i]-d1

                                if d2<=n_doses-1 and d2>=0:
                                        fail[d1,i] = FY[d1,d2]
                                        ratio[d1,i] = log10(Res_array_1[int(d1),int(d2),int(fail[d1,i])]/Res_array_2[int(d1),int(d2),int(fail[d1,i])])


                cmap = plt.get_cmap('jet')

                fig, ax = plt.subplots(1, 1, figsize=(10,7))
                for i in range(len(dose_sum)):
                        ax.plot(ratio[:,i],fail[:,i],label='Dose sum: %r' % round(2*dose_sum[i]/(n_doses-1),2),color=cmap(i/len(dose_sum)))
                        ax.scatter(ratio[:,i],fail[:,i],color=cmap(i/len(dose_sum)))
                
                ax.yaxis.set_major_locator(MaxNLocator(integer=True))
                ax.set_xlabel('Log ratio of resistance frequencies in failure year',fontsize=16)
                ax.set_ylabel('Failure year',fontsize=16)
                ax.axvline(0,linestyle='--',color='k')

                box = ax.get_position()
                ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

                # Put a legend to the right of the current axis
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

                FY_max = np.amax(FY)+1
                FY_min = np.amin(FY[FY>0])
                C2_A = {'Col_mat': FY,'Col_label': 'Failure Year','Col_bds_vec': (0,FY_min,),'alpha': 0.6,'CB': False} # np.concatenate(0,np.linspace(FY_min,FY_max,FY_max-FY_min+1))}

                fig_A = {'xtick': tick_vec, 'ytick': tick_vec}
                Overlay_plotter(Col_attr=C2_A,figure_attributes=fig_A) # Con_attr=Cn_A) # Col_attr=Cl_A)
                ##
                # Cn_A = {'Con_mats': (Z),'Con_inline':(None),'Con_levels': ([0.5])} # 'Con_levels': None,
                Overlay_plotter(Col_attr=C2_A, figure_attributes=fig_A)
                x = np.linspace(0,2,100)
                y = np.zeros((len(dose_sum),len(x)))
                d1_vec = np.nan*np.zeros((len(dose_sum),2*n_doses-2))
                d2_vec = np.nan*np.zeros((len(dose_sum),2*n_doses-2))

                for i in range(len(dose_sum)):
                        y[i,:] = [dose_sum[i]/(n_doses-1) - ii for ii in x]
                        d1_vec[i,:] = [kk/(n_doses-1) for kk in range(2*n_doses-2)]
                        d2_vec[i,:] = [(dose_sum[i] - ii)/(n_doses-1) for ii in range(2*n_doses-2)]

                        # print(d1_vec[i,:],d2_vec[i,:],dose_sum[i]/(n_doses-1))
                        
                        plt.plot(x,y[i,:],label='Dose sum: %r' % round(2*dose_sum[i]/(n_doses-1),2),color=cmap(i/len(dose_sum)))
                        plt.scatter(d1_vec[i,:],d2_vec[i,:],color=cmap(i/len(dose_sum)))
                # plt.legend()
                
               


plt.show()