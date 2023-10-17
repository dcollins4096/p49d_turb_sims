
from GL import *
import simulation


#
# TEB amplitudes, slopes
#
range_ext=dt.extents()
ext_x = [dt.extents() for n in range(3)]
ext_y = [dt.extents() for n in range(3)]
def plot_amps_amps(simlist,amps_or_slopes='amps', suffix='', axis='x'):
    plt.close('all')
    #fig,ax = plt.subplots(1,3, sharex=True,sharey=True,figsize=(12,4))
    x_list = ['density']#,'velocity','Htotal']
    x_list = ['ClTT'+axis]
    y_list = ['ClTT'+axis]#,'ClEE'+axis,'ClBB'+axis]
    y_list = ['ClEE'+axis]#,'ClEE'+axis,'ClBB'+axis]
    #y_list = ['Htotal']
    nplots=1
    fig,ax = plt.subplots(nplots,nplots, figsize=(5.5,5.5))#,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    ax=nar([[ax]])
    axlist=ax.flatten()




    x_agg=[]
    y_agg=[]
    for sim in simlist: #sim_colors.simlist:
        this_sim = simulation.corral[sim]
        this_sim.load()
        simulation.set_colors(this_sim,cmap_name=sim_colors.cmap)

        if amps_or_slopes == 'amps':
            DICT = this_sim.ampsA
            aos = 'amps'
            Aalpha = 'A'
            do_log=True
        else:
            DICT = this_sim.slopesA
            aos = 'slopes'
            Aalpha = "\\alpha"
            do_log=False



        kwargs = {"c":[this_sim.color], "marker":this_sim.marker, "s":this_sim.marker_size}
        for nx,field_x in enumerate(x_list):
            for ny, field_y in enumerate(y_list):
                this_x,this_y=DICT[field_x], DICT[field_y],
                ax[ny][nx].scatter( this_x,this_y, **kwargs)
                ext_x[nx](this_x)
                ext_y[ny](this_y)
                x_agg.append(this_x)
                y_agg.append(this_y)

    if do_log:
        for aaa in axlist:
            aaa.set_yscale('log')
            aaa.set_xscale('log')

    if do_log:
        for aaa in axlist:
            y_minor = matplotlib.ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
            aaa.yaxis.set_minor_locator(y_minor)
            aaa.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    else:
        for aaa in axlist:
            aaa.minorticks_on()

    if 1:
        pfit = np.polyfit(x_agg,y_agg,1)
        ax[0][0].plot( x_agg, pfit[0]*nar(x_agg)+pfit[1])
        print(pfit)
    if 0:
        for n in range(nplots):
            for m in range(nplots):
                ax[m][n].set(xlim=ext_x[n].minmax, ylim=ext_y[m].minmax)
                if m==1:
                    ax[n][m].set(yticks=[])
                if n<2:
                    ax[n][m].set(xticks=[])

    label_dict={'density':r'$%s_\rho$'%Aalpha,
               'velocity':r'$%s_v$'%Aalpha,
                 'Htotal':r'$%s_H$'%Aalpha,
                  'ClTTx':r'$%s_{TT}$'%Aalpha,
                  'ClEEx':r'$%s_{EE}$'%Aalpha,
                  'ClBBx':r'$%s_{BB}$'%Aalpha,
                  'ClTTy':r'$%s_{TT}$'%Aalpha,
                  'ClEEy':r'$%s_{EE}$'%Aalpha,
                  'ClBBy':r'$%s_{BB}$'%Aalpha,
                  'ClTTz':r'$%s_{TT}$'%Aalpha,
                  'ClEEz':r'$%s_{EE}$'%Aalpha,
                  'ClBBz':r'$%s_{BB}$'%Aalpha}

    for n in range(nplots):
        ax[n][-1].yaxis.tick_right()
        ax[n][0].set_ylabel(label_dict[y_list[n]])
        ax[-1][n].set_xlabel(label_dict[x_list[n]])

            

    #fig.tight_layout()
    outname = '%s/prim_vs_TEB_%s%s.pdf'%(dl.plotdir,aos,suffix)
    fig.subplots_adjust(top=0.98)
    fig.savefig(outname)
    print(outname)
