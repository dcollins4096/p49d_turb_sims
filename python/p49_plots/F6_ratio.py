


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
    fig,ax = plt.subplots(1,3, figsize=(12,4), sharex=True,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()

    for nf,field in enumerate(ratios):
        field_top,field_bottom = field
        field_top += LOS
        field_bottom += LOS

        for sim in simlist:
            this_sim=simulation.corral[sim]
            this_sim.load()

            qu = "%s/%s"%(field_top[-3:-1].upper(), field_bottom[-3:-1].upper())
            label = r'$%s~ \hat{%s}$'%(qu, LOS)
            #axlist[nf  ].text(1e-2, 5e-3, label)
            axlist[nf  ].set_title( label)
            this_y = this_sim.avg_spectra[field_top]/this_sim.avg_spectra[field_bottom]
            this_x = this_sim.avg_spectra['k2d']
            axlist[nf  ].plot(this_x, this_y, c=this_sim.color , linestyle=this_sim.linestyle)
            print(nf)
            if nf==2:
                axlist[nf].axhline(0.5,c=[0.5]*3)
                axlist[nf].axhline(1.,c=[0.5]*3)
            
    for a in axlist:
        dt.axbonk(a,xscale='log',yscale=LOG_OR_LIN,xlabel=None,ylabel=None, ylim=[1e-3,1000])
    for a in axlist:
        a.set_xlabel(r'$k/k_{max}$')
    axlist[0].set_ylabel('Ratio')

    outname = '%s/%s_ratio_%s.pdf'%(dl.plotdir,prefix, LOS)
    fig.savefig(outname)
    print(outname)

#
# Ratio amplitudes
#
def plot_amps(simlist, LOS='y'):
    ratio_ext = dt.extents()
    plt.close('all')
    #fig,ax = plt.subplots(1,3, sharex=True,sharey=True,figsize=(12,4))
    fig,ax = plt.subplots(3,2, figsize=(8,4))#,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()

    y_ext=dt.extents()
    for nf,field in enumerate(ratios):
        field_top,field_bottom = field
        field_top += LOS
        field_bottom += LOS
        mean_ratio = 0; n_ratio = 0
        for sim in simlist:
            this_sim=simulation.corral[sim]
            this_sim.load()


            this_Ms = this_sim.Ms_mean
            this_Ma = this_sim.Ma_mean

            kwargs = {"c":this_sim.color, "marker":this_sim.marker}
            
            xvals = this_sim.avg_spectra['k2d']
            fit_range =this_sim.get_fitrange(xvals)
            mask = (xvals > fit_range[0])*(xvals < fit_range[1])
            #this_y = this_sim.avg_spectra[field_top][mask].mean()/this_sim.avg_spectra[field_bottom][mask].mean()
            this_y = this_sim.ampsA[field_top]/this_sim.ampsA[field_bottom]
            y_ext(this_y)
            ax[nf ][0].scatter(this_Ms, this_y,  **kwargs)
            ax[nf ][1].scatter(this_Ma, this_y,  **kwargs)
            ratio_ext(this_y)

            name_top = field_top[-3:-1]
            name_bot = field_bottom[-3:-1]
            ylab = r"$A_{%s}/A_{%s}$"%(name_top,name_bot)
            ax[nf][0].set_ylabel(ylab)
            
            #heres an ugly kludge.
            if sim[0] in '456':
                mean_ratio+=this_y
                n_ratio += 1
        mean_ratio = mean_ratio/n_ratio
        print('Mean Ratio %s %s %f'%(field_top, field_bottom, mean_ratio))
        #if nf == 1:
            #print("EE/BB mean ", mean_ratio)
        #ax[nf][0].plot( [0.45, 2.7], [0.5]*2, c=[0.5]*4)
        ax[nf][0].axhline( 0.5, c=[0.5]*4)
        ax[nf][1].axhline( 0.5, c=[0.5]*4)

    ax[2][0].set_xlabel(r'$M_{\rm{s}}$')
    ax[2][1].set_xlabel(r'$M_{\rm{A}}$')

    if LOG_OR_LIN == 'log':
        for aaa in axlist:
            aaa.set_yscale('log')
    #    if LOS == 'x':
    #        aaa.set_ylim(5e-2,500)
    #    elif LOS in ['y','z']:
    #        pass
    #        #aaa.set_ylim(1e-2,5)
    #        #aaa.set_ylim( ratio_ext.minmax)

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


