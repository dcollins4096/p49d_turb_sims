


from GL import *
from cycler import cycler
import queb3
from queb3 import powerlaw_fit as plfit
import davetools as dt
import h5py
from matplotlib.pyplot import cm
import sim_colors
reload(sim_colors)
reload(queb3)
import get_all_quantities as gaq
reload(gaq)
verbose=False
from collections import defaultdict
all_slopes=defaultdict(list)

#read_stuff reads spectra and average quantities
#Reload to re-read, but don't do that every time.
import read_stuff as rs
spectra_dict = rs.spectra_dict
import read_avg_quan as raq
quan3=raq.quan3


#
# Ratio Spectra
#
myc={}
prefix='straight'
ratios = list(zip(['avg_clee', 'avg_clbb', 'avg_clbb'],['avg_cltt','avg_cltt','avg_clee']))
#ratios = list(zip(['avg_cltt', 'avg_cltt', 'avg_clbb'],['avg_clee','avg_clbb','avg_clee']))
LOG_OR_LIN='log'
#prefix='cross'
#ratios = list(zip(['avg_clte', 'avg_cltb', 'avg_cleb'],['avg_clee','avg_clbb','avg_clee']))
LOSES = 'y'
if 1:
    plt.close('all')
    fig,ax = plt.subplots(1,3, figsize=(12,4), sharex=True,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()
    for LOS in LOSES:

        for nf,field in enumerate(ratios):
            field_top,field_bottom = field

            for sim in sim_colors.simlist:
                qu = "%s/%s"%(field_top[-2:].upper(), field_bottom[-2:].upper())
                label = r'$%s~ \hat{%s}$'%(qu, LOS)
                #axlist[nf  ].text(1e-2, 5e-3, label)
                axlist[nf  ].set_title( label)
                this_y = spectra_dict[LOS][sim].spectra[field_top]/spectra_dict[LOS][sim].spectra[field_bottom]
                axlist[nf  ].plot(spectra_dict[LOS][sim].lcent, this_y, c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim])
                
                #lcent=spectra_dict[LOS][sim].lcent
                #fit_range =spectra_dict[LOS][sim].fit_range
                #mask = (lcent > fit_range[0])*(lcent < fit_range[1])
                #mean_y = spectra_dict[LOS][sim].spectra[field_top][mask].mean()/spectra_dict[LOS][sim].spectra[field_bottom][mask].mean()
                #axlist[nf].plot( fit_range, [mean_y]*2, c=[0.5]*3)
                #if nf == 1:
                #    myc[sim] = mean_y
                #    print(mean_y)

    for a in axlist:
        dt.axbonk(a,xscale='log',yscale=LOG_OR_LIN,xlabel=None,ylabel=None, ylim=[1e-3,1000])
    for a in axlist:
        a.set_xlabel(r'$k/k_{max}$')
    axlist[0].set_ylabel('Ratio')

    outname = '%s/%s_ratio_%s.pdf'%(plotdir,prefix, LOS)
    fig.savefig(outname)
    print(outname)

#
# Ratio amplitudes
#
ratio_ext = dt.extents()
if 1:
    plt.close('all')
    #fig,ax = plt.subplots(1,3, sharex=True,sharey=True,figsize=(12,4))
    fig,ax = plt.subplots(3,2, figsize=(8,4))#,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()

    for LOS in LOSES:
        for nf,field in enumerate(ratios):
            field_top,field_bottom = field
            mean_ratio = 0; n_ratio = 0
            for sim in sim_colors.simlist:


                this_Ms = quan3[sim]['msavg'].mean()
                this_Ma = quan3[sim]['maavg'].mean()

                kwargs = {"c":sim_colors.color[sim], "marker":sim_colors.marker[sim]}
                
                lcent=spectra_dict[LOS][sim].lcent
                fit_range =spectra_dict[LOS][sim].fit_range
                mask = (lcent > fit_range[0])*(lcent < fit_range[1])
                this_y = spectra_dict[LOS][sim].spectra[field_top][mask].mean()/spectra_dict[LOS][sim].spectra[field_bottom][mask].mean()
                mean_ratio+=this_y
                n_ratio += 1
                ax[nf ][0].scatter(this_Ms, this_y,  **kwargs)
                ax[nf ][1].scatter(this_Ma, this_y,  **kwargs)
                ratio_ext(this_y)

                qu = "%s/%s"%(field_top[-2:].upper(), field_bottom[-2:].upper())
                ax[nf][0].set_ylabel(qu)
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
    
        outname = '%s/TEB_%s_%s_ratio_%s.pdf'%(plotdir,prefix, LOS,'measured3')
        fig.savefig(outname)
        print(outname)


