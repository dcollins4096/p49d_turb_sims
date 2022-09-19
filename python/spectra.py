
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
quan3 = rs.quan3

plotdir="%s/PigPen"%os.environ['HOME']
if 1:
    plt.close('all')
    fig,ax = plt.subplots(1,3, figsize=(12,4))#, sharex=True,sharey=True,figsize=(12,4))
    #fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()

    for nf,field in enumerate(['avg_cltt','avg_clee','avg_clbb']):
        for sim in sim_colors.simlist:

            proj=rs.proj_dict['x'][sim]
            label="%s %0.1f"%(sim, spectra_dict['y'][sim].slopes[field])
            axlist[nf].plot(proj.lcent,   spectra_dict['y'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim], label=label)
            axlist[nf].set_title(r'$C_\ell^{%s}$'%( field[-2:].upper()))

            #proj=rs.proj_dict['y'][sim]
            #label="%s %0.1f"%(sim, spectra_dict['y'][sim].slopes[field])
            #axlist[nf+3].plot(proj.lcent, spectra_dict['y'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim], label=label)
            #axlist[nf+3].plot(proj.lcent, spectra_dict['z'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim])
    for a in axlist:
        dt.axbonk(a,xscale='log',yscale='log',xlabel=None,ylabel=None)
#    for a in axlist[:3]:
#        a.legend(loc=0)
    for a in axlist:
        a.set_xlabel(r'$k/k_{max}$')
    for a in [axlist[0]]:#,axlist[3]]:
        a.set_ylabel(r'$P(k)$')

    fig.savefig('%s/alpha_TEB.pdf'%plotdir)

#
# Spectra plots (both in two rows.  Don't use this.)
#

plotdir=os.environ['HOME']+"/PigPen"
if 0:
    plt.close('all')
    fig,ax = plt.subplots(2,3, sharex=True,sharey=True,figsize=(12,4))
    fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()

    for nf,field in enumerate(['avg_cltt','avg_clee','avg_clbb']):
        for sim in sim_colors.simlist:

            proj=rs.proj_dict['x'][sim]
            label="%s %0.1f"%(sim, spectra_dict['x'][sim].slopes[field])
            axlist[nf].plot(proj.lcent,   spectra_dict['x'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim], label=label)

            proj=rs.proj_dict['y'][sim]
            label="%s %0.1f"%(sim, spectra_dict['y'][sim].slopes[field])
            axlist[nf+3].plot(proj.lcent, spectra_dict['y'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim], label=label)
            #axlist[nf+3].plot(proj.lcent, spectra_dict['z'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim])
    for a in axlist:
        dt.axbonk(a,xscale='log',yscale='log',xlabel=None,ylabel=None)
    for a in axlist[:3]:
        a.legend(loc=0)
    for a in axlist[3:]:
        a.set_xlabel(r'$k/k_max$')
    for a in [axlist[0],axlist[3]]:
        a.set_ylabel(r'$P(k)$')

    fig.savefig('%s/alpha_rho_v_H.pdf'%plotdir)


#
# Primitive spectra
#

if 0:
    plt.close('all')
    fig,ax = plt.subplots(1,3, sharex=True,sharey=True,figsize=(12,4))
    fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()

    for nf,field in enumerate(['avg_d','avg_v','avg_h']):
        for sim in sim_colors.simlist:
            label="%s %0.2f"%(sim, spectra_dict['x'][sim].slopes[field])
            print(label)
            kwargs={'c':sim_colors.color[sim], 'linestyle':sim_colors.linestyle[sim], 
                    'label':label}
            axlist[nf].plot(proj.lcent,   spectra_dict['x'][sim].spectra[field], **kwargs)
            kwargs.pop('label')
            pline(spectra_dict['x'][sim], axlist[nf], field, c='k')
            
    for a in axlist:
        dt.axbonk(a,xscale='log',yscale='log',xlabel=None,ylabel=None)
    #for a in axlist[:3]:
    #    a.legend(loc=1)
    for a in axlist[3:]:
        a.set_xlabel(r'$k/k_max$')
    for a in [axlist[0],axlist[2]]:
        a.set_ylabel(r'$P(k)$')

    fig.savefig('%s/spectra_rho_v_H.pdf'%plotdir)

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
