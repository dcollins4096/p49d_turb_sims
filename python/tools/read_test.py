
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
#the import does the reading.
import read_stuff as rs
#reload(read_stuff)  
