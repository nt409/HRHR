



# #----------------------------------------------------------------------------------------------      
def optimal_animator(Res_vec_1,Res_vec_2,Res_vec_1_int,Res_vec_2_int,Yield1,Yield1_int,Yield_array,Y_int,H,H_int,dose_array_success,dose_array_success_int,interp,success_years,n_seasons,animate_int=1,animate_not_int=0):
    label = 'G'
    bounds_vector = (1,)
    if interp==1:
        matrix = H_int[:,:,0]
        numberofdoses_array = len(H_int[:,0,0])
        array= Y_int[:,:,0,0]
        numberofdoses_contour = len(H_int[:,0,0])
    else:
        matrix = H[:,:,0]
        numberofdoses_array = len(H[:,0,0])
        array = (Yield_array[:,:,0,0])
        numberofdoses_contour = len(H[:,0,0])
    labels = ([0,params.Yield_threshold])
    alpha = 1
    title = 'G in year 1'

    fig = plt.figure(figsize=(14,12))
    ax  = fig.add_subplot(221)

    x_lab = 'Fungicide 1 dose'
    y_lab = 'Fungicide 2 dose'

    x = np.linspace(0-0.5/(numberofdoses_array-1),1+0.5/(numberofdoses_array-1),numberofdoses_array+1)
    y = np.linspace(0-0.5/(numberofdoses_array-1),1+0.5/(numberofdoses_array-1),numberofdoses_array+1)
    X, Y = np.meshgrid(x,y)  

    Maximum = ceil(np.amax(matrix))
    bounds    = np.linspace(bounds_vector[0],Maximum,101)
    ticks    = np.linspace(bounds_vector[0],Maximum,Maximum-bounds_vector[0]+1)

    # define the colormap
    cmap, norm = colourmap_fn(bounds)

    ax.pcolormesh(X,Y,np.transpose(matrix),cmap=cmap, norm=norm,alpha = alpha,linewidth=0,antialiased=True) # np.transpose(matrix)
    # antialiased - attempt to remove grid lines that appear due to overlap of alpha=0.5 squares

    ###    
    x1 = np.linspace(0,1,numberofdoses_contour)
    y1 = np.linspace(0,1,numberofdoses_contour)
    X1, Y1 = np.meshgrid(x1,y1)
    ##

    ax.contourf(X1, Y1,np.transpose(array),colors=('k'),alpha = 0.5,levels = labels)#np.transpose(array)
    ax.set_xlim((0,1))
    ax.set_ylim((0,1))
    
    if animate_int==1:
        ax.scatter(2*dose_array_success_int[0,0],2*dose_array_success_int[1,0],color='c')
    if animate_not_int==1:
        ax.scatter(2*dose_array_success[0,0],2*dose_array_success[1,0],color='c')
    
    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal(size="5%", pad=0.6, pack_start=False)
    # create a second axes for the colorbar
    ax2 = fig.add_axes(cax)
    mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=ticks, boundaries=bounds)#, format='%1i')
    ax2.set_ylabel(label)

    ax.set_title(title)
    ax.set_ylabel(y_lab)
    ax.set_xlabel(x_lab)
# #----------------------------------------------------------------------------------------------
    #RF
    ax3 = fig.add_subplot(222)
    values1 = log10(Res_vec_1[0]) # same whether 'int' or not
    values2 = log10(Res_vec_2[0]) # same whether 'int' or not
    values3 = values1+values2

    y_min = min(values1,values2)-1
    y_pos = [1]
    if animate_not_int==1:
        ax3.plot(y_pos,values1,label='log(R_1)',color='b')
        ax3.plot(y_pos,values2,label='log(R_2)',color='r')
        ax3.plot(y_pos,values3,label='log(R_1 R_2)',color='k')
    
    if animate_int==1:
        values1_int = log10(Res_vec_1_int[0])
        values2_int = log10(Res_vec_2_int[0])
        values3_int = values1_int+values2_int

        ax3.plot(y_pos,values1_int,label='Int log(R_1)',color='g')
        ax3.plot(y_pos,values2_int,label='Int log(R_2)',color='y')
        ax3.plot(y_pos,values3_int,label='Int log(R_1 R_2)',color='c')
    
    # ax3.axhline(log10(0.5),linestyle='--',color='k')
    ax3.set_xlim((0,success_years+3))
    ax3.set_ylim((y_min,0))
    ax3.set_ylabel('log_10 resistance freqs, start of year n')
    ax3.set_xlabel('Year')
    ax3.grid()
    ax3.set_title(r'$\Delta$ prod is: ')
# #----------------------------------------------------------------------------------------------
    # Dose
    ax4 = fig.add_subplot(223)
    y_pos = [1]
    if animate_not_int==1:
        values1 = 2*dose_array_success[0,0]
        values2 = 2*dose_array_success[1,0]
        ax4.plot(y_pos,values1,label='D1',color='b')
        ax4.plot(y_pos,values2,label='D2',color='r')
    if animate_int==1:
        values1_int = 2*dose_array_success_int[0,0]
        values2_int = 2*dose_array_success_int[1,0]
        ax4.plot(y_pos,values1_int,label='D1_int',color='g')
        ax4.plot(y_pos,values2_int,label='D2_int',color='y')
    ax4.set_xlim((0,success_years+3))
    ax4.set_ylim((0,1.1))
    ax4.set_ylabel('Doses')
    ax4.set_xlabel('Year')
    ax4.grid()
# #----------------------------------------------------------------------------------------------
    # Yield
    ax5 = fig.add_subplot(224)
    y_min_Y = 90
    y_pos = [1]
    if animate_not_int==1:
        values1 = Yield1[0]
        ax5.plot(y_pos,values1,label='Y1',color='k')
    if animate_int==1:
        values_int = Yield1_int[0]
        ax5.plot(y_pos,values_int,label='Y1_int',color='b')
    ax5.axhline(params.Yield_threshold,linestyle='--',color='k')
    ax5.set_xlim((0,success_years+3))
    ax5.set_ylim((y_min_Y,100))
    ax5.set_ylabel('Yield')
    ax5.set_xlabel('Year')
    ax5.grid()
# #----------------------------------------------------------------------------------------------
    def update(frame_number):
        # Get an index which we can use to re-spawn the oldest raindrop.
        end = min(success_years+1,n_seasons)
        current_index = frame_number % end

        if interp == 1:
            ax.pcolormesh(X,Y,np.transpose(H_int[:,:,current_index]),cmap=cmap, norm=norm,alpha = alpha,linewidth=0,antialiased=True) # np.transpose(matrix)
            # antialiased - attempt to remove grid lines that appear due to overlap of alpha=0.5 squares
            ax.contourf(X1, Y1,np.transpose(Y_int[:,:,0,current_index]),colors=('k'),alpha = 0.7,levels = labels) # np.transpose(matrix)
        else:
            ax.pcolormesh(X,Y,H[:,:,current_index],cmap=cmap, norm=norm,alpha = alpha,linewidth=0,antialiased=True)# np.transpose(matrix)
            # antialiased - attempt to remove grid lines that appear due to overlap of alpha=0.5 squares
            ax.contourf(X1, Y1,np.transpose(Yield_array[:,:,0,current_index]),colors=('k'),alpha = 0.7,levels = labels)# np.transpose(matrix)
        
        if current_index<success_years:    
            if animate_int==1:
                ax.scatter(2*dose_array_success_int[0,current_index],2*dose_array_success_int[1,current_index],color='c')
            if animate_not_int==1:
                ax.scatter(2*dose_array_success[0,current_index],2*dose_array_success[1,current_index],color='c')

        if frame_number<end:
            ax3.clear()
            ax4.clear()
            ax5.clear()
# #----------------------------------------------------------------------------------------------            
            # RF
            y_pos = np.arange(1,current_index+2,1)
            
            if animate_not_int==1:
                values1 = Res_vec_1[0:(current_index+1)]
                values2 = Res_vec_2[0:(current_index+1)]
                values1 = [log10(i) for i in values1]
                values2 = [log10(i) for i in values2]
                values3 = [values1[i]+values2[i] for i in range(len(values1))]

                ax3.plot(y_pos,values1,label='log(R_1)',color='b')
                ax3.plot(y_pos,values2,label='log(R_2)',color='r')
                ax3.plot(y_pos,values3,label='log(R_1 R_2)',color='k')

            if animate_int==1:
                values1_int = Res_vec_1_int[0:(current_index+1)]
                values2_int = Res_vec_2_int[0:(current_index+1)]
                values1_int = [log10(i) for i in values1_int]
                values2_int = [log10(i) for i in values2_int]
                values3_int = [values1_int[i]+values2_int[i] for i in range(len(values1_int))]
    
                ax3.plot(y_pos,values1_int,label='Int log(R_1)',color='g')
                ax3.plot(y_pos,values2_int,label='Int log(R_2)',color='y')
                ax3.plot(y_pos,values3_int,label='Int log(R_1 R_2)',color='c')

            ax3.set_xlim((0,success_years+3))
            ax3.set_ylim((y_min,0))
            # ax3.axhline(log10(0.5),linestyle='--',color='k')
            ax3.set_ylabel('log_10 resistance freqs, start of year n')
            ax3.set_xlabel('Year')
            ax3.legend()
            ax3.grid()
            if len(y_pos)>=2:
                if animate_int==1:
                    diff = exp(values3_int[-1]) - exp(values3_int[-2])
                    diff2 = exp(values3_int[-1])/exp(values3_int[-2])
                if animate_not_int==1:
                    diff = exp(values3[-1]) - exp(values3[-2])
                    diff2 = exp(values3[-1])/exp(values3[-2])
                title_1 = r'$\Delta$ prod is: '+str(np.around(diff,decimals=4))+r', $\Delta$ ratio is: '+str(np.around(diff2,decimals=2))
            else:
                title_1 = r'$\Delta$ prod is: '
            ax3.axvline(success_years,linestyle='--',color='r',alpha=0.5)
            ax3.set_title(title_1)
# #----------------------------------------------------------------------------------------------            
            # dose
            y_pos = np.arange(1,current_index+2,1)
            if animate_not_int==1:
                values1 = 2*dose_array_success[0,0:(current_index+1)]
                values2 = 2*dose_array_success[1,0:(current_index+1)]
                ax4.plot(y_pos,values1,label='D1',color='b')
                ax4.plot(y_pos,values2,label='D2',color='r')
            if animate_int==1:
                values1_int = 2*dose_array_success_int[0,0:(current_index+1)]
                values2_int = 2*dose_array_success_int[1,0:(current_index+1)]
                ax4.plot(y_pos,values1_int,label='D1_int',color='g')
                ax4.plot(y_pos,values2_int,label='D2_int',color='y')
            ax4.set_xlim((0,success_years+3))
            ax4.set_ylim((0,1.1))
            ax4.set_ylabel('Doses')
            ax4.set_xlabel('Year')
            ax4.axvline(success_years,linestyle='--',color='r',alpha=0.5)
            ax4.legend()
            ax4.grid()
# #----------------------------------------------------------------------------------------------
            # Yield
            y_pos = np.arange(1,current_index+2,1)
            if animate_not_int==1:
                values1 = Yield1[0:(current_index+1)]
                ax5.plot(y_pos,values1,label='Y1',color='k')
            if animate_int==1:
                values_int = Yield1_int[0:(current_index+1)]
                ax5.plot(y_pos,values_int,label='Y1_int',color='b')
            ax5.set_xlim((0,success_years+3))
            ax5.set_ylim((y_min_Y,100))
            ax5.axhline(params.Yield_threshold,linestyle='--',color='k')
            ax5.axvline(success_years,linestyle='--',color='r',alpha=0.5)
            ax5.set_ylabel('Yield')
            ax5.set_xlabel('Year')
            ax5.legend()
            ax5.grid()
# #----------------------------------------------------------------------------------------------      
        title = 'G in year %s' % (current_index+1)
        ax.set_title(title)
        return None

    # Construct the animation, using the update function as the animation director.
    FuncAnimation(fig, update, interval=1000)#, blit=False)
    plt.show()
    return None




