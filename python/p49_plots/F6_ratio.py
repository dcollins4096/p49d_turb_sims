


from GL import *

import simulation


#
# Ratio Spectra
#
myc={}
prefix='straight'
ratios = list(zip(['ClEE', 'ClBB', 'ClBB'],['ClTT','ClTT','ClEE']))
LOG_OR_LIN='log'
LOSES = 'y'
def plot_ratios(simlist, LOS='y'):

    plt.close('all')
    fig,axes = plt.subplots(1,3,figsize=(8,4))
    axlist=axes.flatten()

    for nf,field in enumerate(ratios):
        field_top,field_bottom = field
        field_top += LOS
        field_bottom += LOS

        for sim in simlist:
            this_sim=simulation.corral[sim]
            this_sim.load()

            qu = "%s/%s"%(field_top[-3:-1].upper(), field_bottom[-3:-1].upper())
            label = r'$%s~ \hat{%s}$'%(qu, LOS)
            label = r'$%s$'%(qu)
            #axlist[nf  ].text(1e-2, 5e-3, label)
            this_y = this_sim.avg_spectra[field_top]/this_sim.avg_spectra[field_bottom]
            this_x = this_sim.avg_spectra['k2d']
            axlist[nf  ].plot(this_x, this_y, c=this_sim.color , linestyle=this_sim.linestyle)
            if nf==2:
                axlist[nf].axhline(0.5,c=[0.5]*3, linewidth=0.1)
                axlist[nf].axhline(1.,c=[0.5]*3, linewidth=0.1)

            thax=axlist[nf]
            fitrange = this_sim.get_fitrange(this_x)
            thax.axvline(fitrange[0], c=[0.5]*4,linewidth=0.1)
            thax.axvline(fitrange[1], c=[0.5]*4,linewidth=0.1)
            thax.set(ylabel=label)
            
    for a in axlist:
        a.set(xscale='log',yscale=LOG_OR_LIN,xlabel=None, ylim=[1e-2,100])
    for a in axlist:
        a.set_xlabel(r'$k$')

    outname = '%s/%s_ratio_%s.pdf'%(dl.plotdir,prefix, LOS)
    fig.tight_layout()
    fig.savefig(outname)
    print(outname)

#
# Ratio amplitudes
#
def plot_amps(simlist, LOS='y', use_alf=False):
    ratio_ext = dt.extents()
    plt.close('all')
    #fig,ax = plt.subplots(1,3, sharex=True,sharey=True,figsize=(12,4))
    if use_alf:
        ncol=2
        xsize=8
    else:
        ncol=1
        xsize=5
    fig,ax = plt.subplots(3,ncol, figsize=(xsize,4))#,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()

    y_ext=dt.extents()
    for nf,field in enumerate(ratios):
        field_top,field_bottom = field
        field_top += LOS
        field_bottom += LOS
        mean_ratio = 0; n_ratio = 0
        upper_collector=[]
        for sim in simlist:
            this_sim=simulation.corral[sim]
            this_sim.load()


            name_top = field_top[-3:-1]
            name_bot = field_bottom[-3:-1]
            ylab = r"$A_{%s}/A_{%s}$"%(name_top,name_bot)
            this_Ms = this_sim.Ms_mean
            this_Ma = this_sim.Ma_mean

            kwargs = {"c":[this_sim.color], "marker":this_sim.marker, "s":this_sim.marker_size*20}
            
            xvals = this_sim.avg_spectra['k2d']
            fit_range =this_sim.get_fitrange(xvals)
            mask = (xvals > fit_range[0])*(xvals < fit_range[1])
            #this_y = this_sim.avg_spectra[field_top][mask].mean()/this_sim.avg_spectra[field_bottom][mask].mean()
            this_y = this_sim.ampsA[field_top]/this_sim.ampsA[field_bottom]
            y_ext(this_y)
            if use_alf:
                ax[nf ][0].scatter(this_Ms, this_y,  **kwargs)
                ax[nf ][1].scatter(this_Ma, this_y,  **kwargs)
                ax[nf][0].set_ylabel(ylab)
            else:
                ax[nf ].scatter(this_Ms, this_y,  **kwargs)
                ax[nf].set_ylabel(ylab)
            ratio_ext(this_y)

            
            #heres an ugly kludge.
            if sim[0] in '456':
                upper_collector.append(this_y)
                mean_ratio+=this_y
                n_ratio += 1
        print('upper_collector',np.mean(upper_collector), np.std(upper_collector))
        error = np.std(upper_collector)
        mean_ratio = mean_ratio/n_ratio
        print('Mean Ratio %s %s %.2f %.2f'%(field_top, field_bottom, mean_ratio, error))
        #if nf == 1:
            #print("EE/BB mean ", mean_ratio)
        #ax[nf][0].plot( [0.45, 2.7], [0.5]*2, c=[0.5]*4)
        if use_alf:
            ax[nf][0].axhline( 0.5, c=[0.5]*4)
            ax[nf][1].axhline( 0.5, c=[0.5]*4)
            ax[2][0].set_xlabel(r'$M_{\rm{s}}$')
            ax[2][1].set_xlabel(r'$M_{\rm{A}}$')
        else:
            if nf == 2:
                ax[2].axhline( 0.5, c=[0.5]*4)
            ax[2].set_xlabel(r'$M_{\rm{s}}$')


    if LOG_OR_LIN == 'log':
        for aaa in axlist:
            aaa.set_yscale('log')
    #    if LOS == 'x':
    #        aaa.set_ylim(5e-2,500)
    #    elif LOS in ['y','z']:
    #        pass
    #        #aaa.set_ylim(1e-2,5)
    #        #aaa.set_ylim( ratio_ext.minmax)

    if use_alf:
        for n in range(3):
            ax[n][1].set_ylim( ax[n][0].get_ylim())
            ax[n][1].yaxis.tick_right()

    for aaa in axlist:
        y_minor = matplotlib.ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
        aaa.yaxis.set_minor_locator(y_minor)
        aaa.set_ylim(y_ext.minmax)

    outname = '%s/TEB_%s_%s_ratio_%s.pdf'%(dl.plotdir,prefix, LOS,'measured3')
    fig.savefig(outname)
    print(outname)


