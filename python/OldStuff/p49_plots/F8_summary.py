


from GL import *
import queb3
reload(queb3)
import davetools as dt
import sim_colors
reload(sim_colors)

#read_stuff reads spectra and average quantities
#Reload to re-read, but don't do that every time.
import read_stuff as rs
spectra_dict = rs.spectra_dict
import read_avg_quan as raq
quan3=raq.quan3

LOS = 'x'
fig,ax=plt.subplots(1,1)
plotdir =  "/home/dccollins/PigPen"

for sim in sim_colors.simlist:
    the_x=spectra_dict[LOS][sim].amps['avg_clbb']/spectra_dict[LOS][sim].amps['avg_clee']
    the_y=spectra_dict[LOS][sim].slopes['avg_clee']
    kwargs = {"c":sim_colors.color[sim], "marker":sim_colors.marker[sim]}
    ax.scatter( the_x, the_y, **kwargs)
    ax.set(xlabel=r'$A^{BB}/A^{EE}$', ylabel=r'$\alpha_{EE}$')

ax.axhline(-2.45,c= [0.5]*4)
ax.axvline(0.5,c= [0.5]*4)
fig.tight_layout()
outname = '%s/summary.pdf'%plotdir
fig.savefig(outname)
print(outname)
