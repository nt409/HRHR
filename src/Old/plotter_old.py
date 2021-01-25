#-------------------------------------------------------------------
def Contour_plotter(matrix,labels,inline,x_lab=None,y_lab=None,xtick=None,ytick=None,title = None):
    # print(np.ndim(matrix))
    if np.ndim(matrix) == 2: # one matrix to plot
        numberofdoses_x = np.shape(matrix)[0]
        numberofdoses_y = np.shape(matrix)[1]
    if np.ndim(matrix) == 3: # more than one matrix to plot
        numberofdoses_x = np.shape(matrix)[1]
        numberofdoses_y = np.shape(matrix)[2]

    if x_lab is None:
        x_lab = 'Fungicide 1 dose'
    if y_lab is None:
        y_lab='Fungicide 2 dose'
    
    figCS = plt.figure(figsize=(6,6))
    ax  = figCS.add_subplot(111)
    x = np.linspace(0,1,numberofdoses_x)
    y = np.linspace(0,1,numberofdoses_y)
    X, Y = np.meshgrid(x,y)

    
    if np.ndim(matrix) == 2: # one matrix to plot
        if len(labels) >0 and labels is not None:
            CS = plt.contour(X, Y, np.transpose(matrix),levels = labels) # np.transpose(matrix)
        else:
            CS = plt.contour(X, Y, np.transpose(matrix))# np.transpose(matrix)

        if inline=='inline':
            plt.clabel(CS,inline=1,fontsize=10)
        elif inline is None:
            None
        else:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="6%", pad=0.2)
            CB = figCS.colorbar(CS, cax=cax, extend='both')
            cax.set_ylabel(inline)
    if np.ndim(matrix) == 3: # more than one matrix to plot
        for i in range(len(matrix)):
            if len(labels) >0 and labels[i] is not None:
                CS = plt.contour(X, Y, np.transpose(matrix[i]),levels = labels[i])# np.transpose(matrix)
            else:
                CS = plt.contour(X, Y, np.transpose(matrix[i]))# np.transpose(matrix) # no levels here
           
            if inline[i]=='inline':
                plt.clabel(CS,inline=1,fontsize=10)
            elif inline[i] is None:
                None
            else:
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="6%", pad=0.2)
                CB = figCS.colorbar(CS, cax=cax, extend='both')
                cax.set_ylabel(inline[i])
    ax.set_xlim((0,1))
    ax.set_ylim((0,1))
    if ytick is not None:
        ax.set_yticks(np.linspace(0,1,len(ytick)))
        ax.set_yticklabels(ytick)
        # start, end = ax.get_ylim()
        # ax.yaxis.set_ticks(np.arange(start, end, y_dist))
        # ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
    if xtick is not None:
        ax.set_xticks(np.linspace(0,1,len(xtick)))
        ax.set_xticklabels(xtick)
        # start, end = ax.get_xlim()
        # ax.xaxis.set_ticks(np.arange(start, end, x_dist))
        # ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
    ax.set_xlabel(x_lab)
    ax.set_ylabel(y_lab)
    if title is not None:
        ax.set_title(title)
    return None

#-------------------------------------------------------------------
def Colour_plotter(matrix,label,bounds_vector,x_lab=None,y_lab=None,xtick=None,ytick=None,title = None):
    if np.ndim(matrix) == 2:
        # one matrix to plot
        numberofdoses_x = np.shape(matrix)[0]
        numberofdoses_y = np.shape(matrix)[1]
    if np.ndim(matrix) == 3:
        # more than one matrix to plot
        numberofdoses_x = np.shape(matrix)[1]
        numberofdoses_y = np.shape(matrix)[2]

    if x_lab is None:
        x_lab = 'Fungicide 1 dose'
    if y_lab is None:
        y_lab='Fungicide 2 dose'

    fig = plt.figure(figsize=(6,5))
    ax  = fig.add_subplot(111)
    
    x = np.linspace(0-0.5/(numberofdoses_x-1),1+0.5/(numberofdoses_x-1),numberofdoses_x+1)
    y = np.linspace(0-0.5/(numberofdoses_y-1),1+0.5/(numberofdoses_y-1),numberofdoses_y+1)
    X, Y = np.meshgrid(x,y)  

    Maximum = ceil(np.amax(matrix))
    if len(bounds_vector) == 1:
        bounds    = np.linspace(bounds_vector[0],Maximum,101)
        ticks    = np.linspace(bounds_vector[0],Maximum,Maximum-bounds_vector[0]+1)
    if len(bounds_vector) == 2:
        bd1 = [bounds_vector[1]-1]
        bd3    = np.linspace(bounds_vector[1],Maximum,101)
        bounds = np.concatenate((bd1,bd3))
        ticks2  = np.linspace(bounds_vector[1],Maximum,Maximum-bounds_vector[1]+1)
        ticks = np.concatenate(([0],ticks2))

    # define the colormap
    cmap = plt.cm.hot
    # extract all colors from the .hot map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    
    # force the first color entry to be black
    if len(bounds_vector) == 2:
        cmaplist[0] = [0.2, 0.2, 0.2]
    
    # create the new map
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

    # define the bins and normalize
    norm   = mpl.colors.BoundaryNorm(bounds, cmap.N)

    im = ax.pcolormesh(X,Y,np.transpose(matrix),cmap=cmap, norm=norm)# np.transpose(matrix)

    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal(size="5%", pad=0.6, pack_start=False)
    # create a second axes for the colorbar
    ax2 = fig.add_axes(cax)
    cb  = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=ticks, boundaries=bounds)#, format='%1i')
    
    if title is not None:
        ax.set_title(title)
    dose_ticks = [0,0.2,0.4,0.6,0.8,1]
    ax.set_xticks(dose_ticks)
    ax.set_yticks(dose_ticks)
    ax.set_ylabel(y_lab)
    ax.set_xlabel(x_lab)
    if ytick is not None:
        ax.set_yticks(np.linspace(0,1,len(ytick)))
        ax.set_yticklabels(ytick)
    if xtick is not None:
        ax.set_xticks(np.linspace(0,1,len(xtick)))
        ax.set_xticklabels(xtick)
    ax2.set_ylabel(label)
    return None

#----------------------------------------------------------------------------------------------
def plot_Y_RF_and_LY(S_vec,S_vec2,R_vec1_lh,R_vec1_hl,R_vec1_m,R_vec2_lh,R_vec2_hl,R_vec2_m,Y_vec_lh,Y_vec_hl,Y_vec_m,j_lh,j_hl,j_m):
    fig2 = plt.figure(figsize=(16,4))
    ax  = fig2.add_subplot(141)
    ax.axhline(params['Yield_threshold'],linestyle='--',color='k')
    ax.plot(S_vec,Y_vec_m,color='r', label='Mixture')
    ax.plot(S_vec,Y_vec_hl,color='g', label='Alt High/Low')
    ax.plot(S_vec,Y_vec_lh, color='b', label='Alt Low/High')
    #ax.legend(loc='best')
    ax.set_xlabel('Season Number')
    ax.set_ylabel('Yield (% of disease free)')
    ax.grid()
#----------------------------------------------------------------------------------------------
    #fig2 = plt.figure()
    ax2   =fig2.add_subplot(142)
    ax2.plot(S_vec2,R_vec1_m,color='r')#, label='Mixture')
    ax2.plot(S_vec2,R_vec1_hl,color='g')#, label='Alt High/Low')
    ax2.plot(S_vec2,R_vec1_lh, color='b')#, label='Alt Low/High')
    #ax2.legend(loc='best')
    ax2.set_xlabel('Season Number')
    ax2.set_ylabel('Resistance Frequency (Fungicide 1)')
    ax2.grid()
#----------------------------------------------------------------------------------------------
    #fig2 = plt.figure()
    ax3   =fig2.add_subplot(143)
    ax3.plot(S_vec2,R_vec2_m,color='r')#, label='Mixture')
    ax3.plot(S_vec2,R_vec2_hl,color='g')#, label='Alt High/Low')
    ax3.plot(S_vec2,R_vec2_lh, color='b')#, label='Alt Low/High')
    #ax3.legend(loc='best')
    ax3.set_xlabel('Season Number')
    ax3.set_ylabel('Resistance Frequency (Fungicide 2)')
    ax3.grid()
#----------------------------------------------------------------------------------------------
    #fig3 = plt.figure()
    ax4   =fig2.add_subplot(144)
    #ax3.plot(S_vec,R_vec,color='g')#, label='Resistance Frequency')
    y = [j_m,j_hl,j_lh]
    N = len(y)
    x = range(N)
    #ax4.bar(x,y)#F_y_m,color='r')
    N = 3
    ind = np.arange(N)
    width = 0.9
    vals = [j_m,j_hl,j_lh]
    colors = ['r','g','b']
    barlist = ax4.bar(ind, vals, width)
    barlist[0].set_color('r')
    barlist[1].set_color('g')
    barlist[2].set_color('b')
    ax4.set_ylabel('Lifetime yield (multiple of disease-free yield)')
    ax4.grid()
    fig2.legend()#loc='best'
    #plt.show()
    return None

#----------------------------------------------------------------------------------------------
def plot_FYY_SR_LY(Y_v_m,Y_v_lh,Y_v_hl,D_vec,S_v1_m, S_v1_lh, S_v1_hl, S_v2_m, S_v2_lh, S_v2_hl, D_vec2, LY_v_m, LY_v_lh, LY_v_hl, D_vec3):
    fig3 = plt.figure(figsize=(16,4))
    ax  = fig3.add_subplot(141)
    ax.axhline(params['Yield_threshold'],linestyle='--',color='k')
    ax.plot(D_vec,Y_v_m,  color='r', label='Mixture')
    ax.plot(D_vec,Y_v_hl, color='g', label='Alt High/Low')
    ax.plot(D_vec,Y_v_lh, color='b', label='Alt Low/High')
    #ax.legend(loc='best')
    ax.set_xlabel('Dose (Fungicide 2)')
    ax.set_ylabel('First Season Yield (% of disease free)')
    ax.grid()
#----------------------------------------------------------------------------------------------
    #fig2 = plt.figure()
    ax2   =fig3.add_subplot(142)
    ax2.plot(D_vec2,S_v1_m ,color='r')#, label='Mixture')
    ax2.plot(D_vec2,S_v1_hl,color='g')#, label='Alt High/Low')
    ax2.plot(D_vec2,S_v1_lh,color='b')#, label='Alt Low/High')
    #ax2.legend(loc='best')
    ax2.set_xlabel('Dose (Fungicide 2)')
    ax2.set_ylabel('Selection Ratio (Fungicide 1)')
    ax2.grid()
#----------------------------------------------------------------------------------------------
    #fig2 = plt.figure()
    ax3   =fig3.add_subplot(143)
    ax3.plot(D_vec2,S_v2_m ,color='r')#, label='Mixture')
    ax3.plot(D_vec2,S_v2_hl,color='g')#, label='Alt High/Low')
    ax3.plot(D_vec2,S_v2_lh,color='b')#, label='Alt Low/High')
    #ax3.legend(loc='best')
    ax3.set_xlabel('Dose (Fungicide 2)')
    ax3.set_ylabel('Selection Ratio (Fungicide 2)')
    ax3.grid()
#----------------------------------------------------------------------------------------------
    #fig3 = plt.figure()
    ax4   =fig3.add_subplot(144)
    ax4.plot(D_vec3,LY_v_m ,color='r')#, label='Mixture')
    ax4.plot(D_vec3,LY_v_hl,color='g')#, label='Alt High/Low')
    ax4.plot(D_vec3,LY_v_lh, color='b')#, label='Alt Low/High')
    ax4.set_ylabel('Lifetime yield (multiple of disease-free yield)')
    ax4.set_xlabel('High Risk Dose')
    ax4.grid()
    fig3.legend()#loc='best'
    #plt.show()
    return None

#----------------------------------------------------------------------------------------------
def plot_LY_bydose(M,numberofdoses,title,label,bounds,ticks):
    fig4 = plt.figure(figsize=(6,5))
    ax  = fig4.add_subplot(111)

    x = np.linspace(0-0.5/(numberofdoses-1),1+0.5/(numberofdoses-1),numberofdoses+1)
    y = np.linspace(0-0.5/(numberofdoses-1),1+0.5/(numberofdoses-1),numberofdoses+1)
    X, Y = np.meshgrid(x,y)  
    #fig4.colorbar(im, ax=ax)#me
    #fig4.suptitle(title) #was in

    # define the colormap
    cmap = plt.cm.hot
    # extract all colors from the .jet map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    # force the first color entry to be black
    cmaplist[0] = [0.2, 0.2, 0.2]
    # create the new map
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

    # define the bins and normalize
    norm   = mpl.colors.BoundaryNorm(bounds, cmap.N)

    #plot-me
    im = ax.pcolormesh(X,Y,M,cmap=cmap, norm=norm)

    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal(size="5%", pad=0.6, pack_start=False)
    # create a second axes for the colorbar
    ax2 = fig4.add_axes(cax)#[0.9, 0.1, 0.03, 0.8])
    cb  = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=ticks, boundaries=bounds)#, format='%1i')
    #fig.add_axes(cax)
    ax.set_title(title)
    dose_ticks = [0,0.2,0.4,0.6,0.8,1]
    ax.set_xticks(dose_ticks)
    ax.set_yticks(dose_ticks)
    ax.set_ylabel('Fungicide 1 Dose')#swap?
    ax.set_xlabel('Fungicide 2 Dose')#swap?
    ax2.set_ylabel(label)#, size=12)
    return None

#----------------------------------------------------------------------------------------------
def plot_TY_bydose(M,numberofdoses,title,label):
    fig4 = plt.figure(figsize=(7,5))
    ax  = fig4.add_subplot(111)
    
    x = np.linspace(0-0.5/(numberofdoses-1),1+0.5/(numberofdoses-1),numberofdoses+1)
    y = np.linspace(0-0.5/(numberofdoses-1),1+0.5/(numberofdoses-1),numberofdoses+1)
    X, Y = np.meshgrid(x,y)  
    
    #fig4.suptitle(title) #was in

    # define the colormap
    #cmap = plt.cm.hot
    # extract all colors from the .jet map
    #cmaplist = [cmap(i) for i in range(cmap.N)]
    # force the first color entry to be black
    #cmaplist[0] = [0.2, 0.2, 0.2]
    # create the new map
    #cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

    # define the bins and normalize
    #norm   = mpl.colors.BoundaryNorm(bounds, cmap.N)

    #plot-me
    im = ax.pcolormesh(X,Y,M,cmap='hot')#, norm=norm)
    
    #fig4.colorbar(surf, shrink=0.5, aspect=5)
    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal(size="5%", pad=0.6, pack_start=False)
    # create a second axes for the colorbar
    #ax2 = fig4.add_axes(cax)#[0.9, 0.1, 0.03, 0.8])
    ##cb  = mpl.colorbar.ColorbarBase(ax2, cmap='hot', spacing='proportional')#norm=norm, ticks=ticks, boundaries=bounds)#, format='%1i')
    cb = fig4.colorbar(im, ax=ax)#me
    ax.set_title(title)
    dose_ticks = [0,0.2,0.4,0.6,0.8,1]
    ax.set_xticks(dose_ticks)
    ax.set_yticks(dose_ticks)
    ax.set_ylabel('Fungicide 1 Dose')
    ax.set_xlabel('Fungicide 2 Dose')
    cb.set_label(label)#, size=12)
    return None



#----------------------------------------------------------------------------------------------
def plot_LY_bydose_3D(M,numberofdoses,title,label,bounds,ticks):
    x = np.linspace(0-0*0.5/(numberofdoses-1),1+0*0.5/(numberofdoses-1),numberofdoses)
    y = np.linspace(0-0*0.5/(numberofdoses-1),1+0*0.5/(numberofdoses-1),numberofdoses)
    X, Y = np.meshgrid(x,y)  
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_wireframe(X,Y,M)# plot_surface, cmap='hot')# linewidth=0, antialiased=False) # rstride=1, cstride=1, 
    ax.set_title(label)
    return None

#----------------------------------------------------------------------------------------------
def plot_LY_bydose_3D_scatter(M,numberofdoses,title,label,bounds,ticks):
    x = np.linspace(0-0*0.5/(numberofdoses-1),1+0*0.5/(numberofdoses-1),numberofdoses)
    y = np.linspace(0-0*0.5/(numberofdoses-1),1+0*0.5/(numberofdoses-1),numberofdoses)
    X, Y = np.meshgrid(x,y)  
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    surf = ax.scatter(X,Y, M)# plot_surface, cmap='hot')# linewidth=0, antialiased=False) # rstride=1, cstride=1, 
    ax.set_title('Lifetime Yield by Dose')
    return None
    

def plot_doses(dose_11_vec1,dose_12_vec1,dose_21_vec1,dose_22_vec1,dose_11_vec2,dose_12_vec2,dose_21_vec2,dose_22_vec2,dose_11_vec3,dose_12_vec3,dose_21_vec3,dose_22_vec3,dose_11_vec4,dose_12_vec4,dose_21_vec4,dose_22_vec4,dose_11_vec5,dose_12_vec5,dose_21_vec5,dose_22_vec5,S_vec):
    fig = plt.figure(figsize=(18,3))
    ax  = fig.add_subplot(151)
    ax.plot(S_vec,dose_11_vec1,color='r', label='Spray 1,1')
    ax.plot(S_vec,dose_12_vec1,color='y', label='Spray 1,2')
    ax.plot(S_vec,dose_21_vec1,color='g', label='Spray 2,1')
    ax.plot(S_vec,dose_22_vec1,color='b', label='Spray 2,2')
    ax.set_xlabel('Season Number')
    ax.set_ylabel('Dose, system 1')
    ax.set_ylim(0,1)
    ax.grid()
#----------------------------------------------------------------------------------------------
    ax2  = fig.add_subplot(152)
    ax2.plot(S_vec,dose_11_vec2,color='r')#, label='Spray 1,1')
    ax2.plot(S_vec,dose_12_vec2,color='y')#, label='Spray 1,2')
    ax2.plot(S_vec,dose_21_vec2,color='g')#, label='Spray 2,1')
    ax2.plot(S_vec,dose_22_vec2,color='b')#, label='Spray 2,2')
    ax2.set_xlabel('Season Number')
    ax2.set_ylabel('Dose, system 2')
    ax2.set_ylim(0,1)
    ax2.grid()
#----------------------------------------------------------------------------------------------
    ax3  = fig.add_subplot(153)
    ax3.plot(S_vec,dose_11_vec3,color='r')#, label='Spray 1,1')
    ax3.plot(S_vec,dose_12_vec3,color='y')#, label='Spray 1,2')
    ax3.plot(S_vec,dose_21_vec3,color='g')#, label='Spray 2,1')
    ax3.plot(S_vec,dose_22_vec3,color='b')#, label='Spray 2,2')
    ax3.set_xlabel('Season Number')
    ax3.set_ylabel('Dose, system 3')
    ax3.set_ylim(0,1)
    ax3.grid()
#----------------------------------------------------------------------------------------------
    ax4  = fig.add_subplot(154)
    ax4.plot(S_vec,dose_11_vec4,color='r')#, label='Spray 1,1')
    ax4.plot(S_vec,dose_12_vec4,color='y')#, label='Spray 1,2')
    ax4.plot(S_vec,dose_21_vec4,color='g')#, label='Spray 2,1')
    ax4.plot(S_vec,dose_22_vec4,color='b')#, label='Spray 2,2')
    ax4.set_xlabel('Season Number')
    ax4.set_ylabel('Dose, system 4')
    ax4.set_ylim(0,1)
    ax4.grid()
#----------------------------------------------------------------------------------------------
    ax5  = fig.add_subplot(155)
    ax5.plot(S_vec,dose_11_vec5,color='r')#, label='Spray 1,1')
    ax5.plot(S_vec,dose_12_vec5,color='y')#, label='Spray 1,2')
    ax5.plot(S_vec,dose_21_vec5,color='g')#, label='Spray 2,1')
    ax5.plot(S_vec,dose_22_vec5,color='b')#, label='Spray 2,2')
    ax5.set_xlabel('Season Number')
    ax5.set_ylabel('Dose, system 5')
    ax5.set_ylim(0,1)
    ax5.grid()
    fig.subplots_adjust(wspace = 0.4)
    fig.legend()
    return None

#----------------------------------------------------------------------------------------------
def plot_innoc(I_vec1,I_vec2,I_vec3,I_vec4,I_vec5,S_vec):
    fig = plt.figure(figsize=(18,3))
    ax  = fig.add_subplot(151)
    ax.plot(S_vec,I_vec1,color='b')
    ax.set_xlabel('Season Number')
    ax.set_ylabel('Innoculum, system 1')
    #ax.set_ylim(0,1)
    ax.grid()
#----------------------------------------------------------------------------------------------
    ax2  = fig.add_subplot(152)
    ax2.plot(S_vec,I_vec2,color='b')
    ax2.set_xlabel('Season Number')
    ax2.set_ylabel('Innoculum, system 2')
    #ax2.set_ylim(0,1)
    ax2.grid()
#----------------------------------------------------------------------------------------------
    ax3  = fig.add_subplot(153)
    ax3.plot(S_vec,I_vec3,color='b')#, label='Spray 1,1')
    ax3.set_xlabel('Season Number')
    ax3.set_ylabel('Innoculum, system 3')
    #ax3.set_ylim(0,1)
    ax3.grid()
#----------------------------------------------------------------------------------------------
    ax4  = fig.add_subplot(154)
    ax4.plot(S_vec,I_vec4,color='r')
    ax4.set_xlabel('Season Number')
    ax4.set_ylabel('Innoculum, system 4')
    #ax4.set_ylim(0,1)
    ax4.grid()
#----------------------------------------------------------------------------------------------
    ax5  = fig.add_subplot(155)
    ax5.plot(S_vec,I_vec5,color='r')#, label='Spray 1,1')
    ax5.set_xlabel('Season Number')
    ax5.set_ylabel('Innoculum, system 5')
    #ax5.set_ylim(0,1)
    ax5.grid()
    fig.subplots_adjust(wspace = 0.4)
    #fig.legend()
    return None

#----------------------------------------------------------------------------------------------
def plot_innoc2(I_vec1,I_vec2,I_vec3,I_vec4,I_vec5,S_vec):
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.plot(S_vec,I_vec1,color='r')
    ax.plot(S_vec,I_vec2,color='y')
    ax.plot(S_vec,I_vec3,color='g')
    ax.plot(S_vec,I_vec4,color='c')
    ax.plot(S_vec,I_vec5,color='b')
    ax.set_xlabel('Season Number')
    ax.set_ylabel('Innoculum')
    #ax.set_ylim(0,1)
    ax.grid()
    fig.legend()
    return None