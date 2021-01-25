import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import fsolve
from scipy.interpolate import RegularGridInterpolator
from math import ceil, log, floor, log10, exp
from Functions_and_plotting.parameters_HRHR import params
from Functions_and_plotting.functions_HRHR import primary_calculator, master_loop_one_tactic
import warnings
#-------------------------------------------------------------------
def default_key_fn(atts,default_atts):
    for key in default_atts.keys():
        if not(key in atts):
            atts[key] = default_atts[key]
    return atts


#-------------------------------------------------------------------
def colourmap_fn(bounds,grey=False):
    cmap = plt.get_cmap('hot')
    # extract all colors from the .hot map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    if grey:
        # force the first color entry to be grey
        cmaplist[0] = [0.4, 0.4, 0.4] # 0.2 0.2 0.2
    # create the new map
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    norm   = mpl.colors.BoundaryNorm(bounds, cmap.N)
    return cmap, norm


#-------------------------------------------------------------------
# asexual plot functions
def invisible_axes(big_ax):
    big_ax.spines['top'].set_color('none')
    big_ax.spines['bottom'].set_color('none')
    big_ax.spines['left'].set_color('none')
    big_ax.spines['right'].set_color('none')
    big_ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    return big_ax


#-------------------------------------------------------------------
def Overlay_plotter(Col_attr={},Con_attr={},Scat_attr={},figure_attributes={}):
    # defaults
    dose_ticks = [0,0.2,0.4,0.6,0.8,1]
    figure_defaults = {'x_lab': 'Fungicide 1 dose','y_lab': 'Fungicide 2 dose','xtick': False,'ytick': False,'xlim': (0,1),'ylim': (0,1),'grid': False, 'title': None, 'fig': None, 'sub': None, 'ax': None}
    Col_defaults = {'Col_mat':None,'Col_label':None,'Col_bds_vec':None,'alpha':1,'CB':True}
    Con_defaults = {'Con_mats':None,'Con_levels':None,'Con_inline':None,'Contourf':None,'alpha':0.6}
    Scat_defaults = {'x_scat':None,'y_scat':None,'x_scat_k':None,'y_scat_k':None,'scat_colours':None,'scat_leg':None,'markers':None}

    defs = [figure_defaults,Col_defaults,Con_defaults,Scat_defaults]
    attrs = [figure_attributes,Col_attr,Con_attr,Scat_attr]
    for default_list, attribute_dictionary in zip(defs,attrs):
        attribute_dictionary = default_key_fn(attribute_dictionary,default_list)
    

    if figure_attributes['ax'] is not None:
        ax = figure_attributes['ax']
    else:
        if figure_attributes['sub'] is None:
            fig = plt.figure(figsize=(7,6))
            ax  = fig.add_subplot(111)
        else:
            fig = figure_attributes['fig']
            ax  = fig.add_subplot(figure_attributes['sub'])
    


    # Colour plot
    if Col_attr['Col_mat'] is not None:
        numberofdoses_col_x = np.shape(Col_attr['Col_mat'])[0]
        numberofdoses_col_y = np.shape(Col_attr['Col_mat'])[1]
        x = np.linspace(0-0.5/(numberofdoses_col_x-1),1+0.5/(numberofdoses_col_x-1),numberofdoses_col_x+1)
        y = np.linspace(0-0.5/(numberofdoses_col_y-1),1+0.5/(numberofdoses_col_y-1),numberofdoses_col_y+1)
        X, Y = np.meshgrid(x,y)
        ##
        Maximum = ceil(np.amax(Col_attr['Col_mat']))
        Minimum = floor(np.amin(Col_attr['Col_mat']))
        bounds  = np.linspace(Minimum,Maximum,101)
        ticks  = np.linspace(Minimum,Maximum,11)
        if Col_attr['Col_bds_vec'] is not None and len(Col_attr['Col_bds_vec']) == 1:
            bounds    = np.linspace(Col_attr['Col_bds_vec'][0],Maximum,101)
            ticks    = np.linspace(Col_attr['Col_bds_vec'][0],Maximum,Maximum-Col_attr['Col_bds_vec'][0]+1)
        if Col_attr['Col_bds_vec'] is not None and len(Col_attr['Col_bds_vec']) == 2:
            bd1 = [Col_attr['Col_bds_vec'][0]]
            bd3    = np.linspace(Col_attr['Col_bds_vec'][1],Maximum,101)
            bounds = np.concatenate((bd1,bd3))
            ticks2  = np.linspace(Col_attr['Col_bds_vec'][1],Maximum,Maximum-floor(Col_attr['Col_bds_vec'][1])+1)
            ticks = np.concatenate((bd1,ticks2))
        if Col_attr['Col_bds_vec'] is not None and len(Col_attr['Col_bds_vec']) == 3:
            bd1 = [Col_attr['Col_bds_vec'][0]]
            bd3    = np.linspace(Col_attr['Col_bds_vec'][1],Col_attr['Col_bds_vec'][2],101)
            bounds = np.concatenate((bd1,bd3))
            ticks2  = np.linspace(Col_attr['Col_bds_vec'][1],Col_attr['Col_bds_vec'][2],ceil(Col_attr['Col_bds_vec'][2])-floor(Col_attr['Col_bds_vec'][1])+1)
            ticks = np.concatenate((bd1,ticks2))
        if Col_attr['Col_bds_vec'] is not None and len(Col_attr['Col_bds_vec']) > 3:
            bounds = Col_attr['Col_bds_vec']
            ticks = Col_attr['Col_bds_vec']
        
        if Col_attr['Col_bds_vec'] is not None and len(Col_attr['Col_bds_vec']) == 2:
            grey = True
        else:
            grey = False

        cmap, norm = colourmap_fn(bounds,grey)

        if Con_attr['Con_mats'] is None:
            ax.pcolormesh(X,Y,np.transpose(Col_attr['Col_mat']),cmap=cmap, alpha = Col_attr['alpha'], norm=norm)
        else:
            ax.pcolormesh(X,Y,np.transpose(Col_attr['Col_mat']),cmap=cmap, norm=norm,alpha = Con_attr['alpha'],linewidth=0,antialiased=True) # np.transpose(Col_mat)
        # antialiased - attempt to remove grid lines that appear due to overlap of alpha=0.5 squares

        divider = make_axes_locatable(ax)
        cax = divider.new_horizontal(size="5%", pad=0.6, pack_start=False)
        # create a second axes for the colorbar
        alpha_CB=1
        if Col_attr['CB']:
            if Col_attr['alpha'] is not None:
                alpha_CB = Col_attr['alpha']
            ax2 = fig.add_axes(cax)
            mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', alpha=alpha_CB, ticks=ticks, boundaries=bounds)#, format='%1i')
            
            if Col_attr['Col_label'] is not None:
                ax2.set_ylabel(Col_attr['Col_label'],fontsize=16)

    # Contours
    if Con_attr['Con_mats'] is not None:
        numberofdoses_contour_x = np.shape(Con_attr['Con_mats'])[-2]
        numberofdoses_contour_y = np.shape(Con_attr['Con_mats'])[-1]
        ###    
        x1 = np.linspace(0,1,numberofdoses_contour_x)
        y1 = np.linspace(0,1,numberofdoses_contour_y)
        X1, Y1 = np.meshgrid(x1,y1)
        ##
        def contour_setup(s_contf,s_inline,s_con_levels,s_con_mats,index = None):
            # initialise
            contf_colours = None
            con_lev = None
            # single
            if index is None:
                con_lev = s_con_levels
                ##
                Contour_criteria = s_contf is None 
                if not(Contour_criteria) and isinstance(s_contf,str):
                    contf_colours = s_contf
                ##
                Contour_f_criteria = s_contf is not None
            # not single, one of several
            if index is not None:
                if s_con_levels is not None and len(s_con_levels)>0 and s_con_levels[index] is not None:
                    con_lev = s_con_levels[index]
                ##
                Contour_criteria = s_contf is None or s_contf[index] is None 
                if not(Contour_criteria) and isinstance(s_contf[index],str):
                    contf_colours = (s_contf[index])
                ##
                Contour_f_criteria = not(s_contf is None or s_contf[index] is None)

            ####
            # plot
            if Contour_criteria:
                CS = ax.contour(X1, Y1, np.transpose(s_con_mats),levels=(con_lev)) # np.transpose(matrix)
            if Contour_f_criteria:
                CS = ax.contourf(X1, Y1, np.transpose(s_con_mats),colors=contf_colours,levels=con_lev) # np.transpose(matrix)

            # inline or CB?
            if s_inline=='inline':
                ax.clabel(CS,inline=1,fontsize=10)
            elif not(s_inline is None):
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="6%", pad=0.2)
                fig.colorbar(CS, cax=cax, extend='both')
                cax.set_ylabel(s_inline)
            return None
        
        if np.ndim(Con_attr['Con_mats']) == 2: # one matrix to plot
            contour_setup(s_contf=Con_attr['Contourf'],
            s_inline = Con_attr['Con_inline'],
            s_con_levels = Con_attr['Con_levels'],
            s_con_mats=Con_attr['Con_mats'])
        ##   
        if np.ndim(Con_attr['Con_mats']) == 3: # more than one matrix to plot
            for i in range(len(Con_attr['Con_mats'])):
                contour_setup(s_contf=Con_attr['Contourf'],
                s_inline = Con_attr['Con_inline'][i],
                s_con_levels = Con_attr['Con_levels'],
                s_con_mats=Con_attr['Con_mats'][i],
                index=i)
            ###
    
    #
    if figure_attributes['ytick'] is False:
        ax.set_yticks(dose_ticks)
    else:
        ax.set_yticks(np.linspace(0,1,len(figure_attributes['ytick'])))
        ax.set_yticklabels(figure_attributes['ytick'])
    #
    if figure_attributes['xtick'] is False:
        ax.set_xticks(dose_ticks)
    else:
        ax.set_xticks(np.linspace(0,1,len(figure_attributes['xtick'])))
        ax.set_xticklabels(figure_attributes['xtick'])
    #
    ax.set_xlim(figure_attributes['xlim'])
    ax.set_ylim(figure_attributes['ylim'])

    if figure_attributes['title'] is not None:
        ax.set_title(figure_attributes['title'])
    #
    if figure_attributes['grid']:
        ax.grid()

    if figure_attributes['sub'] is None:
        ax.set_xlabel(figure_attributes['x_lab'],fontsize=16)
        ax.set_ylabel(figure_attributes['y_lab'],fontsize=16)

    # for m,c,x,y in zip(plot_i['marker_types'],plot_i['colours'],plot_i['scatter_x_data'],plot_i['scatter_y_data']):
    #     a.scatter(x,y,marker=m,color=c)

    if Scat_attr['x_scat'] is not None:
        if Scat_attr['scat_colours'] is None:
            Scat_attr['scat_colours'] = 'c'
        if Scat_attr['markers'] is None:
            Scat_attr['markers'] = '+'
            ax.scatter(Scat_attr['x_scat'],Scat_attr['y_scat'],color= Scat_attr['scat_colours'],marker=Scat_attr['markers'])
        else:
            for i in range(len(Scat_attr['x_scat'])):
                ax.scatter(Scat_attr['x_scat'][i],Scat_attr['y_scat'][i],color= Scat_attr['scat_colours'],marker=Scat_attr['markers'][i],s=100) # scat_colours[i]
    if Scat_attr['x_scat_k'] is not None:
        ax.scatter(Scat_attr['x_scat_k'],Scat_attr['y_scat_k'],color='k',marker='+')
    if Scat_attr['scat_leg'] is not None:
        for i in range(len(Scat_attr['scat_leg'][0])):
            if Scat_attr['scat_leg'][1][i] in Scat_attr['scat_colours']:
                ax.scatter(-1,-1,color=Scat_attr['scat_leg'][1][i],marker='+',label=Scat_attr['scat_leg'][0][i])
        ax.legend()

    plt.tight_layout()

    return None


#----------------------------------------------------------------------------------------------
def bar_plotter(solutiont, solution,title=None,bar4 = None,bar3=None,sub1=None,sub2=None,sub3=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    norm1 = 1/(solution[0,params.PS_ind] +  solution[0,params.PRS_ind] + solution[0,params.PSR_ind] + solution[0,params.PR_ind])
    norm2 = 1/(solution[-1,params.IS_ind] + solution[-1,params.IRS_ind] + solution[-1,params.ISR_ind] + solution[-1,params.IR_ind])
    if title is not None:
        ax.set_title(title)
    # else:
    #     ax.set_title('Total disease = '+str(np.around(1/norm2,decimals=4)))


    def bar_setup(objects,values1,values2,values3=None,pos=None,bar_width = 0.35,normalise1=norm1,normalise2=norm2,normalise3=1):
        
        y_pos = np.arange(len(objects))

        values1 = [normalise1*i for i in values1]
        values2 = [normalise2*i for i in values2]
        
        if values3 is not None:
            bar_width = 0.2
            values3 = [normalise3*i for i in values3]

        # left, bottom, width, height
        if pos is not None:
            axes = fig.add_axes(pos)
            # plot the zoomed portion
            axes.bar(y_pos,values1,bar_width,align='center',alpha=0.5,color='b') 
            axes.bar(y_pos+bar_width,values2,bar_width,align='center',alpha=0.5,color='r')
            if values3 is not None:
                axes.bar(y_pos+2*bar_width,values3,bar_width,align='center',alpha=0.5,color='g')

        else:
            plt.bar(y_pos,values1,bar_width,align='center',alpha=0.5,color='b',label=r'Primary, Year $n$')
            plt.bar(y_pos+bar_width,values2,bar_width,align='center',alpha=0.5,color='r',label='Final')
            if values3 is not None:
                plt.bar(y_pos+2*bar_width,values3,bar_width,align='center',alpha=0.5,color='g',label=r'Primary, Year $n+1$')
        if values3 is not None:
            plt.xticks(y_pos + bar_width,objects)
        else:        
            plt.xticks(y_pos + 0.5*bar_width,objects)

        return None

    if bar4 is not None:
        bar_setup(objects = ('rr','rs','sr','ss'),
        values1 = [solution[0,params.PR_ind],solution[0,params.PRS_ind],solution[0,params.PSR_ind],solution[0,params.PS_ind]],
        values2 = [solution[-1,params.IR_ind],solution[-1,params.IRS_ind],solution[-1,params.ISR_ind],solution[-1,params.IS_ind]])
        plt.xlabel('Pathogen Strain',fontsize=16)
        plt.ylabel('Frequency',fontsize=16)
        plt.tight_layout()
        
        if sub1 is not None:
            bar_setup(objects = ('rr',),
            values1 = [solution[0,params.PR_ind]],
            values2 = [solution[-1,params.IR_ind]],
            pos = [.3, .2, .3, .3])
            
        if sub2 is not None:
            bar_setup(objects = ('rs','sr'),
            values1 = [solution[0,params.PRS_ind],solution[0,params.PSR_ind]],
            values2 = [solution[-1,params.IRS_ind],solution[-1,params.ISR_ind]],
            pos = [.3, .6, .3, .3])

        if sub3 is not None:
            bar_setup(objects = ('rr','rs','sr'),
            values1 = [solution[0,params.PR_ind],solution[0,params.PRS_ind],solution[0,params.PSR_ind]],
            values2 = [solution[-1,params.IR_ind],solution[-1,params.IRS_ind],solution[-1,params.ISR_ind]],
            pos = [.3, .3, .3, .3])

    if bar3 is not None:
        X_R, X_RS, X_SR, X_S = primary_calculator(res_prop_1=(solution[-1,params.IR_ind]+solution[-1,params.IRS_ind] )*norm2,
        res_prop_2=(solution[-1,params.IR_ind]+solution[-1,params.ISR_ind] )*norm2,
        p_rr=solution[-1,params.IR_ind]*norm2,
        p_rs=solution[-1,params.IRS_ind]*norm2,
        p_sr=solution[-1,params.ISR_ind]*norm2,
        p_ss=solution[-1,params.IS_ind]*norm2)
        # innoc = params.init_den,
        norm3 = 1/(X_R+X_RS+X_SR+X_S)

        bar_setup(objects = ('rr','rs','sr'),
        values1 = [solution[0,params.PR_ind],solution[0,params.PRS_ind],solution[0,params.PSR_ind]],
        values2 = [solution[-1,params.IR_ind],solution[-1,params.IRS_ind],solution[-1,params.ISR_ind]],
        values3 = [X_R,X_RS,X_SR],
        normalise3=norm3)
        plt.xlabel('Pathogen Strain',fontsize=16)
        plt.ylabel('Frequency',fontsize=16)
        
        plt.tight_layout()

        if sub1 is not None:
            bar_setup(objects = ('rr',),
            values1 = [solution[0,params.PR_ind]],
            values2 = [solution[-1,params.IR_ind]],
            values3 = [X_R],
            normalise3=norm3,
            pos = [.35, .65, .1, .18])
    
    fig.legend()
    return None


#----------------------------------------------------------------------------------------------
def line_or_scatter_plotter(plots,figure_attributes={}):
    # defaults
    defaults = {'legend': False, 'label_fontsize': None, 'title_fontsize': None, 'horiz': False, 'share_x':False, 'share_y':False, 'fs': (7,5)}

    figure_attributes = default_key_fn(figure_attributes,defaults)

    if not figure_attributes['horiz']:
        fig, ax = plt.subplots(len(plots), 1, sharex=figure_attributes['share_x'], sharey=figure_attributes['share_y'], figsize=figure_attributes['fs'])
    else:
        fig, ax = plt.subplots(1, len(plots), sharex=figure_attributes['share_x'], sharey=figure_attributes['share_y'], figsize=figure_attributes['fs'])
    
    for a,plot_i in zip(ax,plots):
        # grey plot
        if 'grey_x_box_pos' in plot_i and plot_i['grey_x_box_pos'] is not None and plot_i['grey_y_box_start'] is not None and plot_i['grey_y_box_end'] is not None:
            a.fill_between(plot_i['grey_x_box_pos'],plot_i['grey_y_box_start'],plot_i['grey_y_box_end'], facecolor='grey')

        plot_i_defaults = {'line_y_data': None, 'line_x_data': None,'marker_types': None,'colours': None,'xlab': None,'ylab': None,'xlim':None,'ylim':None}
        plot_i = default_key_fn(plot_i,plot_i_defaults)

        # line plot
        if plot_i['line_y_data'] is not None:
            for j in range(len(plot_i['line_y_data'])):
                y = plot_i['line_y_data'][j]
                # x
                if plot_i['line_x_data'] is None:
                    n = len(y)
                    x = np.linspace(1,n,n)
                else:
                    x = plot_i['line_x_data'][j]
                # colour, lab
                colour, lab = None, None
                if 'line_label' in plot_i and plot_i['line_label'][j] is not None:
                    lab = plot_i['line_label'][j]
                if 'line_colour' in plot_i and plot_i['line_colour'][j] is not None:
                    colour = plot_i['line_colour'][j]
                a.plot(x,y,color = colour,label = lab)
        
        # scatter plot
        if 'scatter_x_data' in plot_i and 'scatter_y_data' in plot_i:
            for m,c,x,y in zip(plot_i['marker_types'],plot_i['colours'],plot_i['scatter_x_data'],plot_i['scatter_y_data']):
                a.scatter(x,y,marker=m,color=c)
        # axis lab
        if plot_i['ylab'] is not None:
            a.set_ylabel(plot_i['ylab'],fontsize = figure_attributes['label_fontsize'])
        if plot_i['xlab'] is not None:
            a.set_xlabel(plot_i['xlab'],fontsize = figure_attributes['label_fontsize'])
        # axis limits
        if plot_i['xlim'] is not None:
            a.set_xlim(plot_i['xlim'])
        if plot_i['ylim'] is not None:
            a.set_ylim(plot_i['ylim'])
        
    
    # L,R,B,T = (0.12,0.9,0.12,0.9)
    # if ('xlab' in figure_attributes and figure_attributes['xlab'] is not None) or ('ylab' in figure_attributes and figure_attributes['ylab'] is not None):
    #     fig.tight_layout(rect=[0, 0.03, 1, 0.95]) #[0.05, 0.05, 0.95, 0.95]) # [0, 0.03, 1, 0.95]
    #     fig.subplots_adjust(left=L,right=R,bottom=B,top=T)
    # if 'xlab' in figure_attributes and figure_attributes['xlab'] is not None:
    #     xlab_ax = fig.add_axes([L, 0.1, R-L-0.06, 0.01]) #([L, 0.08, R-L-0.06, 0.01]) 
    #     xlab_ax = invisible_axes(xlab_ax)
    #     xlab_ax.set_xlabel(figure_attributes['xlab'],fontsize = figure_attributes['label_fontsize'])
    # if 'ylab' in figure_attributes and figure_attributes['ylab'] is not None:
    #     ylab_ax = fig.add_axes([0.1, B, 0.01, T-B]) # [0.08, B, 0.01, T-B]
    #     ylab_ax = invisible_axes(ylab_ax)
    #     ylab_ax.set_ylabel(figure_attributes['ylab'],fontsize = figure_attributes['label_fontsize'])
    
    L,R,B,T = (0.12,0.9,0.12,0.9)
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
    fig.subplots_adjust(left=L,right=R,bottom=B,top=T)


    if ('xlab' in figure_attributes and figure_attributes['xlab'] is not None) or ('ylab' in figure_attributes and figure_attributes['ylab'] is not None):    
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)

    if 'xlab' in figure_attributes and figure_attributes['xlab'] is not None:
        plt.xlabel(figure_attributes['xlab'],fontsize=figure_attributes['label_fontsize'])
    if 'ylab' in figure_attributes and figure_attributes['ylab'] is not None:
        plt.ylabel(figure_attributes['ylab'],fontsize=figure_attributes['label_fontsize'])

    if figure_attributes['legend']:
        if 'leg_pos' in figure_attributes:
            loc = figure_attributes['leg_pos']
        else:
            loc = 'best'
        fig.legend(loc=loc)
    if 'title' in figure_attributes:
        fig.suptitle(figure_attributes['title'],fontsize = figure_attributes['title_fontsize'])
    return None


#----------------------------------------------------------------------------------------------
def mosaic_plot_2d(Array_to_plot,phi_vec_to_plot,phi_rr_vec,inline_vec_to_plot,index=None,fig_type=None,fs = (10,9),title=None,colorbar_on=True,index_vec=range(9)):
    # setup
    L,R,B,T = (0.12,0.9,0.12,0.9)
    # xtick, ytick = np.linspace(min(phi_vec_to_plot),max(phi_vec_to_plot),5), np.linspace(min(phi_vec_to_plot),max(phi_vec_to_plot),5)
    fig, ax0 = plt.subplots(3, 3, sharex=True, sharey=True, figsize=fs)

    if fig_type is None:
        fig_type='col'
        fs = (11,8)
    

    if fig_type=='col':
        # for colour plots
        dist = phi_vec_to_plot[1]-phi_vec_to_plot[0]
        x = np.linspace(phi_vec_to_plot[0]-0.5*dist,phi_vec_to_plot[-1]+0.5*dist,len(phi_vec_to_plot)+1)
    
        if inline_vec_to_plot[-1]>inline_vec_to_plot[1]:
            ticks2  = np.linspace(inline_vec_to_plot[1],inline_vec_to_plot[-1],ceil(inline_vec_to_plot[-1])-floor(inline_vec_to_plot[1])+1)
            ticks = np.concatenate(([inline_vec_to_plot[0]],ticks2))
            bd3    = np.linspace(inline_vec_to_plot[1],inline_vec_to_plot[-1],101)
            bounds = np.concatenate(([inline_vec_to_plot[0]],bd3))
        else:
            ticks  = inline_vec_to_plot
            bounds = np.linspace(inline_vec_to_plot[0],inline_vec_to_plot[-1],101)

        
        if colorbar_on:
            cmap2, norm = colourmap_fn(bounds,grey=True)
        else:
            cmap2, norm = colourmap_fn(bounds,grey=False)

        ##
        # index_maker=False
        # if index_vec is None:
        #     index_vec = range(9)
        #     index_maker=True
        # if index_vec is None:
        #     index_vec = range(9)
        
        for i in range(9):    
            
            # if index_maker:
            #     ii = i % 3
            #     j = floor(i/3)
            # else:
            ii = i % 3
            j = floor(i/3)
            if index_vec[i]<Array_to_plot.shape[0]:
                ax = ax0[j,ii]

                if index is None:
                    matrix = Array_to_plot[index_vec[i],:,:]
                else:
                    matrix = Array_to_plot[index_vec[i],:,:,index]
                im = ax.pcolormesh(x,x,np.transpose(matrix),cmap=cmap2, norm=norm)
                ax.set_title(r'$log_{10} (p_{rr})=$' + '%s' % round(phi_rr_vec[index_vec[i]],3))
                # if j == 2:
                #     ax.set_xlabel('hi x')
                # if ii == 0:
                #     ax.set_ylabel('hi y')

                ax.set_xlim((phi_vec_to_plot[0],phi_vec_to_plot[-1]))
                ax.set_ylim((phi_vec_to_plot[0],phi_vec_to_plot[-1]))
            else:
                ax = ax0[j,ii]
                ax.spines['top'].set_color('none')
                ax.spines['bottom'].set_color('none')
                ax.spines['left'].set_color('none')
                ax.spines['right'].set_color('none')
                ax.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)


            
        ##
        if colorbar_on:
            fig.tight_layout(rect=[0, 0.03, 1, 0.95])
            fig.subplots_adjust(left=L,right=R,bottom=B,top=T)
            cbar_ax = fig.add_axes([R+0.02, B, 0.02, T-B])
            mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap2, norm=norm, spacing='proportional', ticks=ticks, boundaries=bounds) #cb=
    
    
    if fig_type=='con':
        for i in range(9):
            ii = i % 3
            j = floor(i/3)
            if i<Array_to_plot.shape[0]:
                matrix = Array_to_plot[i,:,:]
                x = phi_vec_to_plot
                ax = ax0[j,ii]
                im = ax.contour(x,x,np.transpose(matrix),levels=inline_vec_to_plot)
                ax.clabel(im,inline=1,fontsize=9)
                ax.set_title(r'$p_{rr}=$' + '%s' % round(phi_rr_vec[i],3))
            else:
                ax = ax0[j,ii]
                ax.spines['top'].set_color('none')
                ax.spines['bottom'].set_color('none')
                ax.spines['left'].set_color('none')
                ax.spines['right'].set_color('none')
                ax.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    

        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        fig.subplots_adjust(left=L,right=R,bottom=B,top=T)
    ##

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.xlabel(r'$log_{10} (p_{rs})$',fontsize=16)
    plt.ylabel(r'$log_{10} (p_{sr})$',fontsize=16)

    fig.suptitle(title,fontsize=16)
    return None


#----------------------------------------------------------------------------------------------
def Equal_dose(Yield,Selection_array_1,Res_array_1,season_step,n_doses,n_seasons):
        D = np.linspace(0,1,n_doses)
        SR1 = np.zeros((n_doses,n_seasons))
        Y   = np.zeros((n_doses,n_seasons))
        RF  = np.zeros((n_doses,n_seasons))
        for ii in range(n_doses):
                for jj in range(n_seasons):
                        SR1[ii,jj] = Selection_array_1[ii,ii,jj]
                        Y[ii,jj]   = Yield[ii,ii,jj]
                        RF[ii,jj]   = Res_array_1[ii,ii,jj]
        fig = plt.figure()
        ax = fig.add_subplot(311)
        ax.set_ylabel('Y')
        # did have ax, ax2 labels Dose
        
        ax2 = fig.add_subplot(312)
        ax2.set_ylabel('SR')
        
        ax3 = fig.add_subplot(313)
        ax3.set_ylabel('RF')
        ax3.set_xlabel('Dose')
        # ax.grid()
        # ax2.grid()
        # ax3.grid()
        cmap = plt.get_cmap('jet')
        for jj in range(1,n_seasons,season_step):
                ax.plot(D,Y[:,jj],color = cmap(jj/(n_seasons-1)),label = 'Year %s' % jj )
                ax2.plot(D,SR1[:,jj],color = cmap(jj/(n_seasons-1)))
                ax3.plot(D,RF[:,jj],color = cmap(jj/(n_seasons-1)))
        fig.legend()
        return None


#----------------------------------------------------------------------------------------------
def Season_plotter(colours,labels=None,SRV1=None,SRV2=None,RV1=None,RV2=None,YV=None,RV1_y_lim=None,SRV1_label=None,SRV2_label=None,RV1_label=None):
    fig = plt.figure(figsize=(18,5))
    kk = 0
    j = 1
    for i in (SRV1,SRV2,RV1,RV2,YV):
        if i is not None:
            kk = kk+1
    if YV is not None:
        sub = str(1) + str(kk) + str(j)
        j = j+1
        n_s = max(len(YV[0]),len(YV))
        ax = fig.add_subplot(sub)
        ax.fill_between([-1,n_s*2],-1,params.Yield_threshold, facecolor='grey')
        ax.set_xlim((0,n_s))
        ax.set_ylim((np.amin(YV)-1,100))
        ax.set_xlabel('Season Number',fontsize=16)
        ax.set_ylabel(r'Yield ($\%$ of DF)',fontsize=16)
        # ax.grid()
#----------------------------------------------------------------------------------------------
    if RV1 is not None:
        sub = str(1) + str(kk) + str(j)
        j = j+1
        ax2 = fig.add_subplot(sub)
        ax2.set_xlabel('Season Number',fontsize=16)
        ax2.set_ylabel(r'RF ($F_1$)',fontsize=16)
        # ax2.grid()
        RV1 = np.asarray(RV1)
        n_s = max(len(RV1[0]),len(RV1))
        ax2.set_xlim((0,n_s))
        if RV1_y_lim is not None:
            ax2.set_ylim(RV1_y_lim)
        if RV1_label is not None:
            ax2.set_ylabel(RV1_label)
#----------------------------------------------------------------------------------------------
    if RV2 is not None:
        sub = str(1) + str(kk) + str(j)
        j = j+1
        ax3 = fig.add_subplot(sub)
        ax3.set_xlabel('Season Number',fontsize=16)
        ax3.set_ylabel(r'RF ($F_2$)',fontsize=16)
        n_s = max(len(RV2[0]),len(RV2))
        ax3.set_xlim((0,n_s))
        # ax3.grid()
#----------------------------------------------------------------------------------------------
    if SRV1 is not None:
        sub = str(1) + str(kk) + str(j)
        j = j+1
        ax4 = fig.add_subplot(sub)
        ax4.set_xlabel('Season Number',fontsize=16)
        ax4.set_ylabel(r'SR ($F_1$)',fontsize=16)
        n_s = SRV1.shape[-1] #        n_s = max(len(SRV1[0]),len(SRV1))
        ax4.set_xlim((0,n_s))
        # ax4.grid()
        if SRV1_label is not None:
            ax4.set_ylabel(SRV1_label)
#----------------------------------------------------------------------------------------------
    if SRV2 is not None:
        sub = str(1) + str(kk) + str(j)
        j = j+1
        ax5 = fig.add_subplot(sub)
        ax5.set_xlabel('Season Number',fontsize=16)
        ax5.set_ylabel(r'SR ($F_2$)',fontsize=16)
        n_s = SRV2.shape[-1] # max(len(SRV2[0]),len(SRV2))
        ax5.set_xlim((0,n_s))
        # ax5.grid()
        if SRV2_label is not None:
            ax5.set_ylabel(SRV2_label)
#----------------------------------------------------------------------------------------------
    if len(colours) > 1:
        for i in range(len(colours)):
            if YV is not None:
                S_vec = np.linspace(1,len(YV[i]),len(YV[i]))
                ax.plot(S_vec,YV[i],color=colours[i],label = labels[i])
                fig.legend()
            if RV1 is not None:
                S_vec2 = np.linspace(0,len(RV1[i])-1,len(RV1[i]))
                ax2.plot(S_vec2,RV1[i],color=colours[i])
            if RV2 is not None:
                S_vec2 = np.linspace(0,len(RV2[i])-1,len(RV2[i]))
                ax3.plot(S_vec2,RV2[i],color=colours[i])
            if SRV1 is not None:
                S_vec = np.linspace(1,len(SRV1[i]),len(SRV1[i]))
                ax4.plot(S_vec,SRV1[i],color=colours[i])
            if SRV2 is not None:
                S_vec = np.linspace(1,len(SRV2[i]),len(SRV2[i]))
                ax5.plot(S_vec,SRV2[i],color=colours[i])
    else:
        if YV is not None:
            S_vec = np.linspace(1,len(YV),len(YV))
            ax.plot(S_vec,YV,color=colours,label = labels)
            fig.legend()
        if RV1 is not None:
            S_vec = np.linspace(0,len(RV1)-1,len(RV1))
            ax2.plot(S_vec,RV1, color=colours)
        if RV2 is not None:
            S_vec = np.linspace(0,len(RV2)-1,len(RV2))
            ax3.plot(S_vec,RV2,color=colours)
        if SRV1 is not None:
            S_vec = np.linspace(1,len(SRV1),len(SRV1))
            ax4.plot(S_vec,SRV1,color=colours)
        if SRV2 is not None:
            S_vec = np.linspace(1,len(SRV2),len(SRV2))
            ax5.plot(S_vec,SRV2,color=colours)
    return None


#----------------------------------------------------------------------------------------------
def plot_disease_dynamics(solutiont, solution, title=None,no_R=None):
    
    cmap = plt.get_cmap('jet')
    fig, ax0 = plt.subplots(1, 4, sharex=False, sharey=False, figsize=(14,6))

    ax  = ax0[0]
    ax2 = ax0[1]
    ax3 = ax0[2]
    ax4 = ax0[3]

    #----------------------------------------------------------------------------------------------
    ax.plot(solutiont, solution[:,params.S_ind],color='g',  label=r'$S$')
    ax.plot(solutiont, solution[:,params.R_ind],color='k',  label=r'$R$')
    ax.axvline(params.T_GS32,linestyle='--',color='k')
    ax.axvline(params.T_GS39,linestyle='--',color='k')
    ax.axvline(params.T_GS61,linestyle='--',color='k')
    ax.axvline(params.T_GS87,linestyle='--',color='k')
    ax.text(0.17,1.025,'GS32',horizontalalignment='center',verticalalignment='bottom',fontsize=7,transform=ax.transAxes)#0.19#(T_GS32,0.165,'GS32',horizontalalignment='center',verticalalignment='center',bbox=dict(facecolor='w', alpha=1))
    ax.text(0.32,1.025,'GS39',horizontalalignment='center',verticalalignment='bottom',fontsize=7,transform=ax.transAxes)#0.31#(T_GS39,0.165,'GS39',horizontalalignment='center',verticalalignment='center',bbox=dict(facecolor='w', alpha=1))
    ax.text(0.51,1.025,'GS61',horizontalalignment='center',verticalalignment='bottom',fontsize=7,transform=ax.transAxes)#(T_GS61,0.165,'GS61',horizontalalignment='center',verticalalignment='center',bbox=dict(facecolor='w', alpha=1))
    ax.text(0.95,1.025,'GS87',horizontalalignment='center',verticalalignment='bottom',fontsize=7,transform=ax.transAxes)#(T_GS87,0.165,'GS87',horizontalalignment='center',verticalalignment='center',bbox=dict(facecolor='w', alpha=1))
    ax.set_ylabel('Leaf Area Index (LAI)')

#----------------------------------------------------------------------------------------------
    ax2.axvline(params.T_GS32,linestyle='--',color='k')
    ax2.axvline(params.T_GS39,linestyle='--',color='k')
    ax2.axvline(params.T_GS61,linestyle='--',color='k')
    ax2.axvline(params.T_GS87,linestyle='--',color='k')
    ax2.text(0.17,1.025,'GS32',horizontalalignment='center',verticalalignment='bottom',fontsize=7,transform=ax2.transAxes)#0.19#(T_GS32,0.165,'GS32',horizontalalignment='center',verticalalignment='center',bbox=dict(facecolor='w', alpha=1))
    ax2.text(0.32,1.025,'GS39',horizontalalignment='center',verticalalignment='bottom',fontsize=7,transform=ax2.transAxes)#0.31#(T_GS39,0.165,'GS39',horizontalalignment='center',verticalalignment='center',bbox=dict(facecolor='w', alpha=1))
    ax2.text(0.51,1.025,'GS61',horizontalalignment='center',verticalalignment='bottom',fontsize=7,transform=ax2.transAxes)#(T_GS61,0.165,'GS61',horizontalalignment='center',verticalalignment='center',bbox=dict(facecolor='w', alpha=1))
    ax2.text(0.95,1.025,'GS87',horizontalalignment='center',verticalalignment='bottom',fontsize=7,transform=ax2.transAxes)#(T_GS87,0.165,'GS87',horizontalalignment='center',verticalalignment='center',bbox=dict(facecolor='w', alpha=1))
    
    if no_R is None:
        ax2.plot(solutiont, solution[:,params.ER_ind],color=cmap(0/3),  label=r'$E_{rr}$')
    ax2.plot(solutiont, solution[:,params.ERS_ind],color=cmap(1/3),  label=r'$E_{rs}$')
    ax2.plot(solutiont, solution[:,params.ESR_ind],color=cmap(2/3),  label=r'$E_{sr}$')
    ax2.plot(solutiont, solution[:,params.ES_ind],color=cmap(3/3),  label=r'$E_{ss}$')
    ax2.set_ylabel('LAI')
#----------------------------------------------------------------------------------------------
    ax3.axvline(params.T_GS32,linestyle='--',color='k')
    ax3.axvline(params.T_GS39,linestyle='--',color='k')
    ax3.axvline(params.T_GS61,linestyle='--',color='k')
    ax3.axvline(params.T_GS87,linestyle='--',color='k')
    ax3.text(0.17,1.025,'GS32',horizontalalignment='center',verticalalignment='bottom',fontsize=7,transform=ax3.transAxes)#0.19#(T_GS32,0.165,'GS32',horizontalalignment='center',verticalalignment='center',bbox=dict(facecolor='w', alpha=1))
    ax3.text(0.32,1.025,'GS39',horizontalalignment='center',verticalalignment='bottom',fontsize=7,transform=ax3.transAxes)#0.31#(T_GS39,0.165,'GS39',horizontalalignment='center',verticalalignment='center',bbox=dict(facecolor='w', alpha=1))
    ax3.text(0.51,1.025,'GS61',horizontalalignment='center',verticalalignment='bottom',fontsize=7,transform=ax3.transAxes)#(T_GS61,0.165,'GS61',horizontalalignment='center',verticalalignment='center',bbox=dict(facecolor='w', alpha=1))
    ax3.text(0.95,1.025,'GS87',horizontalalignment='center',verticalalignment='bottom',fontsize=7,transform=ax3.transAxes)#(T_GS87,0.165,'GS87',horizontalalignment='center',verticalalignment='center',bbox=dict(facecolor='w', alpha=1))
    if no_R is None:
        ax3.plot(solutiont, solution[:,params.IR_ind],color=cmap(0/3),  label=r'$I_{rr}$')
    ax3.plot(solutiont, solution[:,params.IRS_ind],color=cmap(1/3),  label=r'$I_{rs}$')
    ax3.plot(solutiont, solution[:,params.ISR_ind],color=cmap(2/3),  label=r'$I_{sr}$')
    ax3.plot(solutiont, solution[:,params.IS_ind],color=cmap(3/3),  label=r'$I_{ss}$')

    ax3.set_ylabel('LAI')
#----------------------------------------------------------------------------------------------
    ax4.plot(solutiont, solution[:,params.Fung1_ind],color='b',  label=r'$C_1$')
    ax4.plot(solutiont, solution[:,params.Fung2_ind],color='g',  label=r'$C_2$')

    ax4.axvline(params.T_GS32,linestyle='--',color='k')
    ax4.axvline(params.T_GS39,linestyle='--',color='k')
    ax4.axvline(params.T_GS61,linestyle='--',color='k')
    ax4.axvline(params.T_GS87,linestyle='--',color='k')
    ax4.text(0.17,1.025,'GS32',horizontalalignment='center',verticalalignment='bottom',fontsize=7,transform=ax4.transAxes)
    ax4.text(0.32,1.025,'GS39',horizontalalignment='center',verticalalignment='bottom',fontsize=7,transform=ax4.transAxes)
    ax4.text(0.51,1.025,'GS61',horizontalalignment='center',verticalalignment='bottom',fontsize=7,transform=ax4.transAxes)
    ax4.text(0.95,1.025,'GS87',horizontalalignment='center',verticalalignment='bottom',fontsize=7,transform=ax4.transAxes)
    ax4.set_ylabel('Concentration')

    # ax.grid()
    # ax2.grid()
    # ax3.grid()
    # ax4.grid()
#----------------------------------------------------------------------------------------------
    if title is not None:
        fig.suptitle(title)
    fig.legend()

    L,R,B,T = (0.12,0.9,0.12,0.9)
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
    fig.subplots_adjust(left=L,right=R,bottom=B,top=T)

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.xlabel('Time (degree days)',fontsize=16)
#----------------------------------------------------------------------------------------------
    return None
#----------------------------------------------------------------------------------------------








# from Optimal_HRHR_asexual_plots
#----------------------------------------------------------------------------------------------
def cube_scatter_plot(Array_to_plot,phi_vec_to_plot,phi_rr_vec):
    # setup
    marker_index = ['o','<','>','v','^','s','d','x','+','.',',','o','<','>','v','^','s','d','x','+','.',',']
    alpha = 0.6
    cmap = plt.get_cmap('jet')

    fig = plt.figure(figsize=(6,8))
    ax = fig.add_subplot(111, projection='3d')
    ##
    for i in range(len(phi_rr_vec)):
        for j in range(len(phi_vec_to_plot)):
            for k in range(len(phi_vec_to_plot)):
                col = cmap(Array_to_plot[i,j,k]/(np.amax(Array_to_plot)))
                marker = marker_index[int(Array_to_plot[i,j,k])]
                ax.scatter(phi_vec_to_plot[j], phi_vec_to_plot[k], phi_rr_vec[i], color = col, s = 50, marker=marker, alpha = alpha) # note have reordered so that rr on z axis

    for ii in range(ceil(np.amax(Array_to_plot))+1):
        ax.scatter(10,10,10, color = cmap(ii/(np.amax(Array_to_plot))), marker = 'o', s = 50, label = 'E.L.=%s' % ii, alpha = alpha) #marker_index[ii]

    # note have reordered so that rr on z axis
    ax.set_zlabel(r'$log_{10}(p_{rr})$',fontsize=16)
    ax.set_xlabel(r'$log_{10}(p_{rs})$',fontsize=16)
    ax.set_ylabel(r'$log_{10}(p_{sr})$',fontsize=16)

    # ax.set_zlabel('Double Resistant',fontsize=16)
    # ax.set_xlabel('Single Resistant 1',fontsize=16)
    # ax.set_ylabel('Single Resistant 2',fontsize=16)

    ax.set_xlim((min(phi_vec_to_plot),max(phi_vec_to_plot)))
    ax.set_ylim((min(phi_vec_to_plot),max(phi_vec_to_plot)))
    ax.set_zlim((min(phi_rr_vec),max(phi_rr_vec)))

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    return None
#----------------------------------------------------------------------------------------------
def surface_plot_3d(Z,phi_vec_to_plot,phi_rr_vec,fs = (6,8),multiplot=True):
    # setup
    c_vec_2 = range(1,min(Z.shape[-1],9))
    cmap = plt.get_cmap('jet')
    fig = plt.figure()
    X,Y = np.meshgrid(phi_vec_to_plot,phi_vec_to_plot)
    
    if multiplot:
        for j in c_vec_2:
            sub = str(33) + str(j) # str(kk) + str(kk)  + str(j) 
            ax = fig.add_subplot(sub,projection='3d')
            ax.plot_surface(X,Y,Z[:,:,j],color = cmap(j/max(c_vec_2)))
            ax.scatter(0,0,10,color = cmap(j/max(c_vec_2)),alpha = 0.5, label = 'F = %s' % j)
            ax.set_xlim((min(phi_vec_to_plot),max(phi_vec_to_plot)))
            ax.set_ylim((min(phi_vec_to_plot),max(phi_vec_to_plot)))
            ax.set_zlim((min(phi_rr_vec),max(phi_rr_vec)))

        fig.legend()
        fig.suptitle('F = contours')
    if not multiplot:
        ax = fig.gca(projection='3d')
        for contour in c_vec_2:
            ax.plot_surface(X,Y,Z[:,:,contour],color = cmap(contour/max(c_vec_2)),alpha = 0.5)
            ax.scatter(0,0,10,color = cmap(contour/max(c_vec_2)),alpha = 0.5, label = 'F = %s' % contour)

        ax.set_xlim((min(phi_vec_to_plot),max(phi_vec_to_plot)))
        ax.set_ylim((min(phi_vec_to_plot),max(phi_vec_to_plot)))
        ax.set_zlim((min(phi_rr_vec),max(phi_rr_vec)))
        ax.legend()
        ax.set_title('F contours')
    return None
#----------------------------------------------------------------------------------------------
def generate_array(F_array_here,phi_vec_here,phi_rr_vec):
    # setup
    x = phi_rr_vec
    y,z = [phi_vec_here for i in range(2)]
    Int_F   = RegularGridInterpolator((x,y,z),F_array_here,bounds_error=False,fill_value=0)
    Z = np.nan*np.ones((len(phi_vec_here),len(phi_vec_here),int(round(np.amax(F_array_here)+1))))
    
    for contour in range(1,int(round(np.amax(F_array_here)+1))):
        for i in range(len(phi_vec_here)-1):
            for j in range(len(phi_vec_here)-1):
                print(min(F_array_here[:,i,j]),max(F_array_here[:,i,j]),contour)
                if min(F_array_here[:,i,j])<=contour and max(F_array_here[:,i,j])>=contour:
                    p1, p2 = phi_vec_here[i], phi_vec_here[j]
                    def equation(y):
                        return float(Int_F([y,p1,p2])) - contour
                    Z[i,j,contour] = fsolve(equation,-0.001)
                    if Z[i,j,contour]>= -0.001 or Z[i,j,contour]<-10:
                        Z[i,j,contour] = np.nan
                        warnings.warn('Contour Warning')
                        print(contour,i,j,'=contour, i, j')
    return Z
#----------------------------------------------------------------------------------------------
def dose_choice_plotter(phi_space_array,doses,Yield,phi_vec_used,phi_rr_vec,title=None,xlab=None,ylab=None,fs=(11,9)):
    if title is None:
        'Dose effect'
    if xlab is None:
        xlab = r'Dose F1'
    if ylab is None:
        ylab = r'Dose F2'
    
    n_y = phi_space_array.shape[1]
    
    fig = plt.figure(figsize=fs)
    ax = fig.add_subplot(221,projection='3d')
    markers = [r'$0$',r'$1$',r'$2$',r'$3$',r'$4$',r'$5$',r'$6$',r'$7$',r'$8$',r'$9$',r'$10$']
    for jj in range(n_y):
        ax.scatter(phi_space_array[1,jj],phi_space_array[2,jj],phi_space_array[0,jj],marker = markers[jj])

    # note have reordered so that rr on z axis
    ax.set_zlabel(r'$log_{10}(p_{rr})$')
    ax.set_xlabel(r'$log_{10}(p_{rs})$')
    ax.set_ylabel(r'$log_{10}(p_{sr})$')
    ax.set_xlim((min(phi_vec_used),max(phi_vec_used)))
    ax.set_ylim((min(phi_vec_used),max(phi_vec_used)))
    ax.set_zlim((min(phi_rr_vec),max(phi_rr_vec)))

    ax2 = fig.add_subplot(222)
    ax2.grid()
    ax2.plot([0,0.5],[0,0.5],alpha=0.5,linestyle=':',color='k')
    for jj in range(n_y-1):
        ax2.scatter(doses[0,jj],doses[1,jj],marker = markers[jj])
    ax2.set_xlim((0,0.5))
    ax2.set_ylim((0,0.5))
    ax2.set_xlabel(xlab)
    ax2.set_ylabel(ylab)

    ax3 = fig.add_subplot(223)
    ax3.grid()
    ax3.axhline(params.Yield_threshold,linestyle=':',color='k')
    ax3.plot(range(len(Yield)),Yield)
    ax3.scatter(range(len(Yield)),Yield)
    ax3.set_xlabel('Year')
    ax3.set_ylabel('Yield')
    ax3.set_ylim((params.Yield_threshold-1,100))

    # Cumulative_Y = np.zeros(len(Yield))
    # Cumulative_Y[0] = Yield[0]
    # for ii in range(len(Cumulative_Y)-1):
    #     Cumulative_Y[ii+1] = Cumulative_Y[ii]+Yield[ii+1]
    ax4 = fig.add_subplot(224)
    ax4.grid()
    ax4.plot(range(doses.shape[1]),doses[0,:],label='D1')
    ax4.plot(range(doses.shape[1]),doses[1,:],label='D2')
    ax4.scatter(range(doses.shape[1]),doses[0,:])
    ax4.scatter(range(doses.shape[1]),doses[1,:])
    ax4.set_xlabel('Year')
    ax4.set_ylabel('Doses')
    ax4.set_ylim((0,0.5))
    ax4.legend()

    fig.suptitle('Dose effect',fontsize=16)
    return None
#----------------------------------------------------------------------------------------------
def phi_space_navigator(X,Dose_int_1,Dose_int_2,log_phi_minimum = -5):
    # log_phi_minimum means global.JSON(log_phi_min)
    X[0] = max(X[0],log_phi_minimum)
    X[1] = max(X[1],log_phi_minimum)
    X[2] = max(X[2],log_phi_minimum)
    dose_1 = Dose_int_1(X)
    dose_2 = Dose_int_2(X)
    PRR = 10**(X[0])
    PRS = 10**(X[1])
    PSR = 10**(X[2])
    PSS = 1-PRR-PRS-PSR
    if dose_1>0 and dose_2>0:
        output = master_loop_one_tactic(dose_1,dose_2,dose_1,dose_2,p_rr=PRR,p_rs=PRS,p_sr=PSR,p_ss=PSS)
    else:
        warnings.warn('Warning, doses are negative')
        print(dose_1[0],dose_2[0])
        output = master_loop_one_tactic([0.5],[0.5],[0.5],[0.5],p_rr=PRR,p_rs=PRS,p_sr=PSR,p_ss=PSS)
    X2 = np.zeros(3)
    X2[0] = log10(output['PRR'][1])

    X2[1] = log10(output['PRS'][1])
    X2[2] = log10(output['PSR'][1])
    doses = [dose_1,dose_2]
    Yield = output['Yield_vec']
    return X2, doses, Yield
#----------------------------------------------------------------------------------------------
