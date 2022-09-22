"""
Make E plot and slope using 2 methods.

*  queb3.simulation_package points to a simulation.  
BoxSize is in units of 64 zones.
*  

"""
from GL import *
from cycler import cycler
import spectra_tools as st
import queb3
reload(queb3)
from queb3 import powerlaw_fit as plfit
import davetools as dt
import h5py
from matplotlib.pyplot import cm
import sim_colors
reload(queb3)
verbose=False
from collections import defaultdict
all_slopes=defaultdict(list)
plt.close('all')
shortprefix="time"
bad_frames=defaultdict(list)
simlist=sim_colors.simlist
framelist=[range(11,31),range(11,31),range(11,31),range(11,31),range(11,31),range(11,31),range(65,85),range(11,31),range(11,31),range(72,91),range(56,75),range(20,40)]
fig,axes=plt.subplots(1,2)
ax0 = axes[0]
ax1 = axes[1]
slopes={}
amps={}
for i in range(12):
    frames=framelist[i]
    simdes=simlist[i]
    sim_dir = "/data/cb1/Projects/P49_EE_BB/%s"%simdes
    plot_dir = "./plots"
    gen_dir = "./plots"
    avg_clee=0
    axes='x'
    longprefix='%s_%s_%s'%(simdes,shortprefix,axes)
    pack = queb3.simulation_package( directory=sim_dir,frames=frames,prefix=longprefix)
    nplots=0


    hspec = 0
    all_slopes=[]
    all_amps=[]
    mean_spectrum = []
    for frame in frames: 
        if 0:
            oober = st.short_oober(sim_dir, frame=frame)
            if len(glob.glob(oober.get_ds_name(frame) )) == 0:
               bad_frames[simdes].append(frame)
               print("BAAAAAAA", bad_frames)
               continue
            try:
                ds=oober.load(frame)
            except:
                bad_frames[simdes].append(frame)
                raise
            print(ds)
        if 0:
            pack.make_htotal_spectra(frame)
        if 1:
            #proj=pack.read_queb(frame=frames[0],ax=axes,bin_style='dx1')
            fname="%s/DD%04d.products/power_Htotal.h5"%(sim_dir,frame)
            if len(glob.glob(fname)) == 0:
                continue

            hspec=dt.dpy( fname , ['k','power'])

            #ax0.plot(hspec[0], hspec[1].real, linestyle=sim_colors.linestyle[simdes], c=sim_colors.color[simdes])
            slopes[simdes] = queb3.slope_package()
            slopes[simdes].ingest("BT", *queb3.powerlaw_fit(hspec[0], hspec[1].real, [1e-2, 8e-2]))
            #ax1.scatter( slopes[simdes].slope["BT"].real, slopes[simdes].amp["BT"].real, color=[0.5]*4,
            #            marker=sim_colors.marker[simdes])
            mean_spectrum.append( hspec[1].real)
            all_slopes.append( slopes[simdes].slope["BT"].real )
            all_amps.append(slopes[simdes].amp["BT"].real)
        ax1.scatter( all_slopes, all_amps,  color=[0.5]*4, #sim_colors.color[simdes],
                    marker=sim_colors.marker[simdes])
    mean_spectrum = nar(mean_spectrum).mean(axis=0)
    ax0.plot(hspec[0], mean_spectrum, linestyle=sim_colors.linestyle[simdes], c=sim_colors.color[simdes])
    ax1.scatter( nar(all_slopes).mean(), nar(all_amps).mean(),c=  sim_colors.color[simdes],
                marker=sim_colors.marker[simdes])

        #ax1.scatter( mean_slope, mean_amp, color=sim_colors.color[simdes],
        #            marker=sim_colors.marker[simdes])

dt.axbonk(ax0,xscale='log',yscale='log', xlabel='k', ylabel=r'$P(||B||)$')
dt.axbonk(ax1,xscale='linear',yscale='log',xlabel=r'$Spectral\ slope\ \alpha$', ylabel=r'$Amplitude$')
fig.savefig('/home/dccollins/PigPen/spectra.pdf')
plt.close(fig)

            #pack.read_htotspectra(frame)
            #fitrange=proj.determine_fit_range()  #something better

