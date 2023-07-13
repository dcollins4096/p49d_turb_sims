
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
#reload(read_stuff)  
import refit3
#spectra_dict = rs.spectra_dict
spectra_dict = refit3.spectra_dict
quan3 = rs.quan3

plot_dir =  "/home/dccollins/PigPen"
plotdir =  "/home/dccollins/PigPen"

if 0:
    nominal_Ma =nar([sim_colors.Ma[sim] for sim in sub])
    nominal_Ms =nar([sim_colors.Ms[sim] for sim in sub])
    ax[0][3].scatter(nominal_Ma.mean(), fits_v[0], c=colors[0])
    AAA['v'].append(fits_v[0])
    MMMa['v'].append(nominal_Ma.mean())
    MMMa['v']=nar(MMMa['v'])
    AAA['v']=nar(AAA['v'])
    pfit = np.polyfit(MMMa['v'],AAA['v'],1)
    ax[0][3].plot(MMMa['v'],pfit[0]*MMMa['v']+pfit[1],c='k')

#
# TEB amplitudes, slopes
#
LOS = 'y'
range_ext=dt.extents()
simlist = ['half_half']#,'1_half','2_half','3_half']#['1_half']#['half_2', 'half_half']#['half_half']#['1_1']#['2_2'] #['half_half'] # ['half_half']#,'2_2']
simlist = sim_colors.simlist
if 1:
    plt.close('all')
    #fig,ax = plt.subplots(1,3, sharex=True,sharey=True,figsize=(12,4))
    fig,ax = plt.subplots(3,2, figsize=(8,4))#,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()

    for sim in simlist: #sim_colors.simlist:

        do_prim=False; do_TEB=False
        do_amp=False; do_slope=False

        if 1:
            F1 = 'avg_cltt'
            F2 = 'avg_clee'
            F3 = 'avg_clbb'
            S1 = 'TT'
            S2 = 'EE'
            S3 = 'BB'
            which_quan = "TEB_%s"%LOS
            do_TEB=True
        else:
            F1 = 'avg_d'
            F2 = 'avg_v'
            F3 = 'avg_h'
            S1 = r'\rho'
            S2 = 'v'
            S3 = 'h'
            which_quan = "prim"
            do_prim=True

        #
        # slopes or amplitudes
        #
        do_log=False
        if 0:
            #amplitudes
            do_log=True
            do_amp=True
            the_y_T = spectra_dict[LOS][sim].amps[F1]
            the_y_E = spectra_dict[LOS][sim].amps[F2]
            the_y_B = spectra_dict[LOS][sim].amps[F3]
            aos = 'amps'
            Aalpha = 'A'
        elif 1:
            #slopes
            do_slope=True
            suffix='AlphaFixed_25'
            suffix='AlphaPeak'
            the_y_T = spectra_dict[LOS][sim].slopes3[F1].AlphaPeak
            the_y_E = spectra_dict[LOS][sim].slopes3[F2].AlphaPeak
            the_y_B = spectra_dict[LOS][sim].slopes3[F3].AlphaPeak
            #the_y_T = spectra_dict[LOS][sim].slope_bucket.avg_slope[F1]
            #the_y_E = spectra_dict[LOS][sim].slope_bucket.avg_slope[F2]
            #the_y_B = spectra_dict[LOS][sim].slope_bucket.avg_slope[F3]
            #the_std_T = spectra_dict[LOS][sim].slope_bucket.std_slope[F1]
            #the_std_E = spectra_dict[LOS][sim].slope_bucket.std_slope[F2]
            #the_std_B = spectra_dict[LOS][sim].slope_bucket.std_slope[F3]
            #the_y_T = spectra_dict[LOS][sim].slopes[F1]
            #the_y_E = spectra_dict[LOS][sim].slopes[F2]
            #the_y_B = spectra_dict[LOS][sim].slopes[F3]
            aos = 'slopes'
            Aalpha = "\\alpha"

        this_Ms = quan3[sim]['msavg'].mean()
        this_Ma = quan3[sim]['maavg'].mean()
        kwargs = {"c":sim_colors.color[sim], "marker":sim_colors.marker[sim]}
        print(sim, sim_colors.color[sim], the_y_T)
        ax[0][0].scatter(this_Ms, the_y_T,  **kwargs)
        print('b', this_Ms)
        ax[1][0].scatter(this_Ms, the_y_E,  **kwargs)
        ax[2][0].scatter(this_Ms, the_y_B,  **kwargs)

        ax[0][1].scatter(this_Ma, the_y_T,  **kwargs)
        ax[1][1].scatter(this_Ma, the_y_E,  **kwargs)
        ax[2][1].scatter(this_Ma, the_y_B,  **kwargs)
        range_ext(the_y_T)
        range_ext(the_y_E)
        range_ext(the_y_B)

    if do_log:
        for aaa in axlist:
            aaa.set_yscale('log')
    ax[0][0].set_ylabel(r'$%s_{\rm{%s}}$'%(Aalpha,S1))
    ax[1][0].set_ylabel(r'$%s_{\rm{%s}}$'%(Aalpha,S2))
    ax[2][0].set_ylabel(r'$%s_{\rm{%s}}$'%(Aalpha,S3))
    ax[2][0].set_xlabel(r'$M_{\rm{s}}$')
    ax[2][1].set_xlabel(r'$M_{\rm{A}}$')

    if do_TEB:
        ax[0][0].set_title("%s %s"%(aos,LOS))


    if (do_prim and do_amp):
        for aaa in axlist:
            aaa.set_ylim([4e-6,4e-2])
    if (do_TEB and do_amp):
        for aaa in axlist:
            aaa.set_ylim([2e-7,2e-1])
    if do_prim and do_slope:
        for aaa in axlist:
            pass
            #aaa.set_ylim([-1.9,-0.9])
    if do_TEB and do_slope:
        for aaa in axlist:
            pass
            #aaa.set_ylim([-5,-1])

    if do_log:
        for aaa in axlist:
            y_minor = matplotlib.ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
            aaa.yaxis.set_minor_locator(y_minor)
            aaa.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    else:
        for aaa in axlist:
            aaa.minorticks_on()

    for n in range(3):
        ax[n][1].yaxis.tick_right()

            

    outname = '%s/%s_mach_%s_%s.pdf'%(plotdir, which_quan,aos,suffix)
    fig.savefig(outname)
    print(outname)
