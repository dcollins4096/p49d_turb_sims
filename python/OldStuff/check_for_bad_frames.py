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
for i in range(12):
    frames=framelist[i]
    simdes=simlist[i]
    sim_dir = "/data/cb1/Projects/P49_EE_BB/%s"%simdes

    for frame in frames: 
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

