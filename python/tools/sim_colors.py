from GL import *
import davetools as dt


mach_label=r'$\mathcal{M}_{\rm{S}}$'
alf_mach_label=r'$\mathcal{M}_{\rm{A}}$'

mach_norm = mpl.colors.Normalize(vmin=0,vmax=7)

cmap=None
cbar=None

def set_cmap(cmap_name):
    global cmap, cbar
    cmap=cmap_name
    cbar = matplotlib.cm.ScalarMappable(mach_norm,cmap=cmap_name)

if cmap is None:
    #cmap='jet'
    set_cmap( 'gist_rainbow')


