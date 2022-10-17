
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

plotdir =  "/home/dccollins/PigPen"


#
# TEB amplitudes, slopes
#
range_ext=dt.extents()
simlist = ['half_half']#,'1_half','2_half','3_half']#['1_half']#['half_2', 'half_half']#['half_half']#['1_1']#['2_2'] #['half_half'] # ['half_half']#,'2_2']
simlist = sim_colors.simlist
ext_x = [dt.extents() for n in range(3)]
ext_y = [dt.extents() for n in range(3)]
def plot_amps_amps(amps_or_slopes='amps', suffix='', axis='x'):
    plt.close('all')
    #fig,ax = plt.subplots(1,3, sharex=True,sharey=True,figsize=(12,4))
    fig,ax = plt.subplots(3,3, figsize=(8,4))#,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()


    x_list = ['avg_cltt','avg_clee','avg_clbb']
    y_list = ['avg_d','avg_v','avg_h']
    nplots=3

    for sim in simlist: #sim_colors.simlist:

        if amps_or_slopes == 'amps':
            DICT = spectra_dict[axis][sim].amps
            aos = 'amps'
            Aalpha = 'A'
            do_log=True
        else:
            DICT = spectra_dict[axis][sim].slopes
            aos = 'slopes'
            Aalpha = "\\alpha"
            do_log=False



        kwargs = {"c":sim_colors.color[sim], "marker":sim_colors.marker[sim]}
        for nx,field_x in enumerate(x_list):
            for ny, field_y in enumerate(y_list):
                ax[ny][nx].scatter( DICT[field_x], DICT[field_y], **kwargs)
                ext_x[nx](DICT[field_x])
                ext_y[ny](DICT[field_y])

    if do_log:
        for aaa in axlist:
            aaa.set_yscale('log')


#   if (do_prim and do_amp):
#       for aaa in axlist:
#           aaa.set_ylim([4e-6,4e-1])
#   if (do_TEB and do_amp):
#       for aaa in axlist:
#           aaa.set_ylim([2e-7,2e-1])
#   if do_prim and do_slope:
#       for aaa in axlist:
#           pass
#           #aaa.set_ylim([-1.9,-0.9])
#   if do_TEB and do_slope:
#       for aaa in axlist:
#           pass
#           #aaa.set_ylim([-5,-1])

    if do_log:
        for aaa in axlist:
            y_minor = matplotlib.ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
            aaa.yaxis.set_minor_locator(y_minor)
            aaa.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    else:
        for aaa in axlist:
            aaa.minorticks_on()

    for n in range(3):
        for m in range(3):
            ax[m][n].set(xlim=ext_x[n].minmax, ylim=ext_y[m].minmax)
            if m==1:
                ax[n][m].set(yticks=[])
            if n<2:
                ax[n][m].set(xticks=[])

    for n in range(nplots):
        ax[n][2].yaxis.tick_right()
        ax[n][0].set_ylabel(y_list[n])
        ax[2][n].set_xlabel(x_list[n])

            

    #fig.tight_layout()
    outname = '%s/amp_amp_%s_%s%s.pdf'%(plotdir, which_quan,aos,suffix)
    fig.savefig(outname)
    print(outname)
plot_amps_amps(amps_or_slopes='slopes',axis='y')
