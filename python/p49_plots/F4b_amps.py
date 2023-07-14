
from GL import *
import queb3
reload(queb3)
import davetools as dt
import sim_colors
reload(sim_colors)
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

plotdir =  "/home/dccollins/PigPen"


#
# TEB amplitudes, slopes
#
range_ext=dt.extents()
simlist = ['half_half']#,'1_half','2_half','3_half']#['1_half']#['half_2', 'half_half']#['half_half']#['1_1']#['2_2'] #['half_half'] # ['half_half']#,'2_2']
simlist = sim_colors.simlist
def plot_amps_slopes(amps_or_slopes='amps',prim_or_teb='prim', suffix='', axis='x'):
    plt.close('all')
    #fig,ax = plt.subplots(1,3, sharex=True,sharey=True,figsize=(12,4))
    fig,ax = plt.subplots(3,2, figsize=(5.5,5.5))#,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()

    for sim in simlist: #sim_colors.simlist:

        do_prim=False; do_TEB=False
        do_amp=False; do_slope=False

        if prim_or_teb=='teb':
            F1 = 'avg_cltt'
            F2 = 'avg_clee'
            F3 = 'avg_clbb'
            S1 = 'TT'
            S2 = 'EE'
            S3 = 'BB'
            which_quan = "TEB_%s"%axis
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
        if amps_or_slopes=='amps':
            #amplitudes
            do_log=True
            do_amp=True
            the_y_T = spectra_dict[axis][sim].amps[F1]
            the_y_E = spectra_dict[axis][sim].amps[F2]
            the_y_B = spectra_dict[axis][sim].amps[F3]
            aos = 'amps'
            Aalpha = 'A'
        elif 1:
            #slopes
            do_slope=True
            the_y_T = spectra_dict[axis][sim].slopes[F1]
            the_y_E = spectra_dict[axis][sim].slopes[F2]
            the_y_B = spectra_dict[axis][sim].slopes[F3]
            aos = 'slopes'
            Aalpha = "\\alpha"

        this_Ms = quan3[sim]['msavg'].mean()
        this_Ma = quan3[sim]['maavg'].mean()
        kwargs = {"c":sim_colors.color[sim], "marker":sim_colors.marker[sim]}
        ax[0][0].scatter(this_Ms, the_y_T,  **kwargs)
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

    #if do_TEB:
    #    ax[0][0].set_title("%s %s"%(aos,axis))


    if (do_prim and do_amp):
        for aaa in axlist:
            aaa.set_ylim([4e-6,4e-1])
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

            

    fig.subplots_adjust(top=0.92)
    #fig.tight_layout()
    if amps_or_slopes=='amps':
        fig.suptitle('Amplitudes', y=0.96)
    else:
        fig.suptitle('Slopes', y=0.96)
    outname = '%s/%s_%s%s.pdf'%(plotdir, which_quan,aos,suffix)
    fig.savefig(outname)
    print(outname)
plot_amps_slopes(amps_or_slopes='amps', prim_or_teb='prim',axis='y')
plot_amps_slopes(amps_or_slopes='slopes', prim_or_teb='prim',axis='y')
plot_amps_slopes(amps_or_slopes='amps',   prim_or_teb='teb',axis='y')
plot_amps_slopes(amps_or_slopes='slopes', prim_or_teb='teb',axis='y')
