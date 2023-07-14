
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
spectra_dict = rs.spectra_dict

plotdir="%s/PigPen"%os.environ['HOME']
def plot_slopes(prim_or_teb='teb',axis='y'):
    if prim_or_teb=='teb':
        product_list = ['avg_cltt','avg_clee','avg_clbb']
        suffix='TEB'
    else:
        product_list = ['avg_d','avg_v','avg_h']
        suffix='Prim'
    plt.close('all')
    fig,ax = plt.subplots(1,3, figsize=(8,3))#, sharex=True,sharey=True,figsize=(12,4))
    #fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()

    proj_to_spit_out=None
    for nf,field in enumerate(product_list):
        for sim in sim_colors.simlist:

            proj=rs.proj_dict[axis][sim]
            proj_to_spit_out=proj
            label="%s %0.1f"%(sim, spectra_dict['y'][sim].slopes[field])
            lcent_to_use = proj.lcent/proj.lcent.max()
            axlist[nf].plot(lcent_to_use,   spectra_dict['y'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim], label=label)
            lab= field[-2:].upper()
            if lab[0]=='_':
                char = {'D':r'\rho','V':'v'}.get(lab[-1],lab[-1])
                lab=char+char
            axlist[nf].set_title(r'$C_\ell^{%s}$'%(lab))

    for a in axlist:
        dt.axbonk(a,xscale='log',yscale='log',xlabel=None,ylabel=None)
#    for a in axlist[:3]:
#        a.legend(loc=0)
    for a in axlist:
        a.set_xlabel(r'$k/k_{max}$')
    for a in [axlist[0]]:#,axlist[3]]:
        a.set_ylabel(r'$P(k)$')

    fig.tight_layout()
    outname='%s/spectra_%s_%s.pdf'%(plotdir,suffix,axis)
    fig.savefig(outname)
    print(outname)
    return proj_to_spit_out

#plot_slopes(prim_or_teb='teb','x')
#plot_slopes(prim_or_teb='prim','x')
proj=plot_slopes(prim_or_teb='teb',axis='y')
plot_slopes(prim_or_teb='prim',axis='y')

#
# Primitive Slope vs TEB slope
#

if 0:
    plt.close('all')
    fig,ax = plt.subplots(3,3, sharex=True,sharey=True,figsize=(12,4))
    #fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()

    for nfx,fieldx in enumerate(['avg_d','avg_v','avg_h']):
        for nfy,fieldy in enumerate(['avg_cltt','avg_clee','avg_clbb']):
            for sim in sim_colors.simlist:
                #label="%s %0.2f"%(sim, spectra_dict['x'][sim].slopes[field])
                #print(label)
                #kwargs={'c':sim_colors.color[sim], 'linestyle':sim_colors.linestyle[sim], 
                #        'label':label}
                ax[nfy][nfx].scatter( spectra_dict['y'][sim].slopes[fieldx],
                                    spectra_dict['y'][sim].slopes[fieldy],
                                   c=sim_colors.color[sim], marker=sim_colors.marker[sim])
                #kwargs.pop('label')
                #pline(spectra_dict['x'][sim], axlist[nf], field, c='k')

            
    for a in axlist:
        dt.axbonk(a,xscale='linear',yscale='linear',xlabel=None,ylabel=None)
    ax[0][0].set_ylabel(r'$\alpha_\rho$')
    ax[1][0].set_ylabel(r'$\alpha_v$')
    ax[2][0].set_ylabel(r'$\alpha_H$')
    ax[2][0].set_xlabel(r'$\alpha_{\rm{TT}}$')
    ax[2][1].set_xlabel(r'$\alpha_{\rm{EE}}$')
    ax[2][2].set_xlabel(r'$\alpha_{\rm{BB}}$')


    fig.savefig('%s/slope_TTEEBB_vs_rho_v_H.pdf'%plotdir)
