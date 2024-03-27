from GL import *
import davetools as dt


mach_label=r'$\mathcal{M}_{\rm{S}}$'
alf_mach_label=r'$\mathcal{M}_{\rm{A}}$'

mach_norm = mpl.colors.Normalize(vmin=0,vmax=7)

cmap=None
cbar=None

planck_E_slope = -2.42 #LR71
planck_B_slope = -2.54 #LR71
planck_ratio = 0.53 #0.48 to 0.53 for patches.  0.5 is LR42, 0.53 is LR71
planck_TE = 0.355
planck_TB = 0.05


def set_cmap(cmap_name):
    global cmap, cbar
    cmap=cmap_name
    cbar = matplotlib.cm.ScalarMappable(mach_norm,cmap=cmap_name)

if cmap is None:
    #cmap='jet'
    set_cmap( 'gist_rainbow')


