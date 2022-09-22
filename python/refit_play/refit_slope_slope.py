
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
    fig,ax = plt.subplots(3,2, figsize=(6,8))#,sharey=True)
    axlist=ax.flatten()
    for aaa in axlist:
        aaa.set_aspect('equal')


    e_t = dt.extents()
    e_e = dt.extents()
    e_b = dt.extents()
    e_d = dt.extents()
    e_v = dt.extents()
    e_h = dt.extents()
    for sim in simlist: #sim_colors.simlist:

        print('====== ',sim)
        do_prim=False; do_TEB=False
        do_amp=False; do_slope=False

        for ppp in [0,1]:
            if ppp==0:
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
            the_y_T = spectra_dict[LOS][sim].slopes3[F1].AlphaPeak
            the_y_E = spectra_dict[LOS][sim].slopes3[F2].AlphaPeak
            the_y_B = spectra_dict[LOS][sim].slopes3[F3].AlphaPeak
            the_x_T = spectra_dict[LOS][sim].slopes[F1]
            the_x_E = spectra_dict[LOS][sim].slopes[F2]
            the_x_B = spectra_dict[LOS][sim].slopes[F3]
            aos = 'slopes'
            Aalpha = "\\alpha"
            print("S1",S1)
            ylab1=r'$%s_{\rm{%s}}$'%(Aalpha,S1)
            print(ppp,ylab1)
            ax[0][ppp].set_ylabel( ylab1)
            ax[1][ppp].set_ylabel(r'$%s_{\rm{%s}}$'%(Aalpha,S2))
            ax[2][ppp].set_ylabel(r'$%s_{\rm{%s}}$'%(Aalpha,S3))
            ax[0][ppp].scatter( the_x_T, the_y_T,marker=sim_colors.marker[sim],color=sim_colors.color[sim])
            ax[1][ppp].scatter( the_x_E, the_y_E,marker=sim_colors.marker[sim],color=sim_colors.color[sim])
            ax[2][ppp].scatter( the_x_B, the_y_B,marker=sim_colors.marker[sim],color=sim_colors.color[sim])
            if ppp==0:
                e_t(nar([ the_x_T, the_y_T]))
                e_e(nar([ the_x_E, the_y_E]))
                e_b(nar([ the_x_B, the_y_B]))
            else:
                e_d(nar([ the_x_T, the_y_T]))
                e_v(nar([ the_x_E, the_y_E]))
                e_h(nar([ the_x_B, the_y_B]))



    ax[0][0].plot( e_t.minmax, e_t.minmax, c=[0.5]*4)
    ax[1][0].plot( e_e.minmax, e_e.minmax, c=[0.5]*4)
    ax[2][0].plot( e_b.minmax, e_b.minmax, c=[0.5]*4)
    ax[0][1].plot( e_d.minmax, e_d.minmax, c=[0.5]*4)
    ax[1][1].plot( e_v.minmax, e_v.minmax, c=[0.5]*4)
    ax[2][1].plot( e_h.minmax, e_h.minmax, c=[0.5]*4)
    if do_log:
        for aaa in axlist:
            aaa.set_yscale('log')

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

            

    outname = '%s/refit_slope_slope.png'%plotdir
    print(outname)
    fig.savefig(outname)
