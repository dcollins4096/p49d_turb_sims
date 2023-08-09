from GL import *
import davetools as dt


mach_norm = mpl.colors.Normalize(vmin=0,vmax=7)
if 'cmap' not in dir():
    #cmap='jet'
    cmap = 'gist_rainbow'
cbar = matplotlib.cm.ScalarMappable(mach_norm,cmap=cmap)
mach_label=r'$\mathcal{M}_{\rm{S}}$'
alf_mach_label=r'$\mathcal{M}_{\rm{A}}$'

