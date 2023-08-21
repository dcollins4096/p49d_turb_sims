
from GL import *

import simulation
import simulation_info.all_sims

def plot_avg_spectra(simlist,prim_or_teb='teb',axis='y'):

    if prim_or_teb == 'prim':
        product_list = ['density','velocity','Htotal']
        axis_label=""
    else:
        product_list = ['ClTT'+axis,'ClEE'+axis,'ClBB'+axis]
        axis_label="_%s"%axis

    fig,axes = plt.subplots(1,3,figsize=(8,4))

    for ns,sim in enumerate(simlist):
        this_sim = simulation.corral[sim]
        this_sim.load()
        for np,prod in enumerate(product_list):
            thax=axes[np]
            if prim_or_teb == 'prim':
                xvals = this_sim.avg_spectra['k3d']
            else:
                xvals = this_sim.avg_spectra['k2d']
            thax.plot( xvals, this_sim.avg_spectra[prod],c=this_sim.color,linestyle=this_sim.linestyle)
            title_dict={'density':r'$C_\ell^{\rho\rho}$', 'velocity':r'$C_\ell^{vv}$', 'Htotal':r'$C_\ell^{HH}$'}
            title_dict['ClTTy']=r'$C_\ell^{TT}\hat{y}$'
            title_dict['ClBBy']=r'$C_\ell^{BB}\hat{y}$'
            title_dict['ClEEy']=r'$C_\ell^{EE}\hat{y}$'
            thax.set(xscale='log',yscale='log',xlabel='k',title=title_dict[prod])

            fitrange = this_sim.get_fitrange(xvals)
            #print(nar(fitrange)/xvals[1])
            #print(xvals/xvals[1])
            thax.axvline(fitrange[0], c=[0.5]*4,linewidth=0.1)
            thax.axvline(fitrange[1], c=[0.5]*4,linewidth=0.1)


    outname='%s/multi_spectra_%s%s.pdf'%(dl.plotdir,prim_or_teb,axis_label)
    fig.tight_layout()
    fig.savefig(outname)
    print(outname)

"""
old stuff.

def old_plot_slopes(simlist,prim_or_teb='teb',axis='y'):
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
#proj=plot_slopes(prim_or_teb='teb',axis='y')
#plot_slopes(prim_or_teb='prim',axis='y')

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
    """
