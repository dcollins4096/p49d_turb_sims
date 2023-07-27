

from GL import *

import simulation


#
# reb, ttb, rte spectra
#
def plot_hist2(simlist,LOS='y'):
    plt.close('all')

    ext=dt.extents()
    #for nf,field in enumerate(['r_TE'+LOS,'r_TB'+LOS,'r_EB'+LOS]):
    #    for sim in simlist:
    #        this_sim=simulation.corral[sim]
    #        this_sim.load()
    #        q= np.abs(this_sim.avg_spectra[field])
    #        ext(q[q>0])
    fig,ax=plt.subplots(1,1)

    for sim in simlist:
        this_sim=simulation.corral[sim]
        this_sim.load()
        the_x=this_sim.avg_spectra['r_TB'+LOS]
        the_y=this_sim.avg_spectra['r_EB'+LOS]
        ax.plot(the_x,the_y,c=this_sim.color,linestyle=this_sim.linestyle)



    outname = '%s/h2_XY.pdf'%dl.plotdir
    fig.savefig(outname)
    print(outname)

#
def plot_meanvar(simlist,LOS='y'):
    plt.close('all')
    fig,ax = plt.subplots(1,3, sharex=True,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()

    for nf,field in enumerate(['r_TE'+LOS,'r_TB'+LOS,'r_EB'+LOS]):
        for sim in simlist:
            this_sim=simulation.corral[sim]
            this_sim.load()
            qu = field[-3:-1].upper()
            label = r'$r_{%s}$'%qu
            axlist[nf].set_title(field)
            xvals = this_sim.avg_spectra['k2d']
            fit_range =this_sim.get_fitrange(xvals)
            mask = (xvals > fit_range[0])*(xvals < fit_range[1])
            signal = this_sim.avg_spectra[field][mask]
            mean = signal.mean()
            std  = signal.std()
            if mean < 0:
                ec='k'
            else:
                ec=None
            axlist[nf].scatter(np.abs(mean),std, c=this_sim.color,marker=this_sim.marker,edgecolor=ec)
            axlist[nf].set(xlabel=r'$\mu$', ylabel=r'$\sigma$')
            if nf==0:
                axlist[nf].axvline(0.35,c=[0.5]*3)
            
            #axlist[nf+3].plot(spectra_dict['y'][sim].lcent, spectra_dict['y'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim])
            #axlist[nf+6].plot(proj.lcent, spectra_dict['z'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim])
    for a in axlist:
        #dt.axbonk(a,xscale='log',yscale='log',xlabel=None,ylabel=None)
        #a.set_yscale('symlog',linthresh=0.09)
        a.set_yscale('linear')
        #a.set_ylim([-1,1])
        #a.set_ylim([-10,10])
        a.set_xscale('symlog',linthresh=1e-2)
#    for a in axlist:
#        a.set_xlabel(r'$r_{XY}$')
#        a.axhline( 5e-2,c=[0.5]*4)
#        a.axhline(-5e-2,c=[0.5]*4)
#    for a in [axlist[0]]:
#        a.set_ylabel(r'$r_{XY}$')
#        a.set_yticks([-1,-0.1,-0.01,0.01,0.1,1])

    outname = '%s/musigma_XY.pdf'%dl.plotdir
    fig.savefig(outname)
    print(outname)


def plot_hist(simlist,LOS='y'):
    plt.close('all')
    fig,ax = plt.subplots(1,3, sharex=True,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()

    ext=dt.extents()
    for nf,field in enumerate(['r_TE'+LOS,'r_TB'+LOS,'r_EB'+LOS]):
        for sim in simlist:
            this_sim=simulation.corral[sim]
            this_sim.load()
            q= np.abs(this_sim.avg_spectra[field])
            ext(q[q>0])

    sym_bins=False
    if sym_bins:
        bins1 = np.geomspace(ext.minmax[0],ext.minmax[1],64)
        print(bins1)
        bins2 = np.concatenate([-bins1[::-1],bins1])
    else:
        bins2=np.geomspace(ext.minmax[0],ext.minmax[1],64)
    for nf,field in enumerate(['r_TE'+LOS,'r_TB'+LOS,'r_EB'+LOS]):
        for sim in simlist:
            this_sim=simulation.corral[sim]
            this_sim.load()
            qu = field[-3:-1].upper()
            label = r'$r_{%s} \parallel$'%qu
            #axlist[nf  ].text(0.1, -0.8, label)
            label = r'$r_{%s} \perp$'%qu
            axlist[nf].set_title(field)
            #axlist[nf+3].text(0.1, -0.8, label)
            xvals = this_sim.avg_spectra['k2d']
            fit_range =this_sim.get_fitrange(xvals)
            mask = (xvals > fit_range[0])*(xvals < fit_range[1])
            the_h = this_sim.avg_spectra[field]
            if sym_bins == False:
                the_h = np.abs(this_sim.avg_spectra[field])
            axlist[nf].hist(the_h,histtype='step',color=this_sim.color,linestyle=this_sim.linestyle,bins=bins2)
            axlist[nf].set(xlabel=field)
            
            #axlist[nf+3].plot(spectra_dict['y'][sim].lcent, spectra_dict['y'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim])
            #axlist[nf+6].plot(proj.lcent, spectra_dict['z'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim])
    for a in axlist:
        #dt.axbonk(a,xscale='log',yscale='log',xlabel=None,ylabel=None)
        #a.set_yscale('symlog',linthresh=0.09)
        a.set_yscale('linear')
        #a.set_ylim([-1,1])
        #a.set_ylim([-10,10])
        a.set_xscale('symlog',linthresh=1e-2)
#    for a in axlist:
#        a.set_xlabel(r'$r_{XY}$')
#        a.axhline( 5e-2,c=[0.5]*4)
#        a.axhline(-5e-2,c=[0.5]*4)
#    for a in [axlist[0]]:
#        a.set_ylabel(r'$r_{XY}$')
#        a.set_yticks([-1,-0.1,-0.01,0.01,0.1,1])

    outname = '%s/h_XY.pdf'%dl.plotdir
    fig.savefig(outname)
    print(outname)

def plot_spectra(simlist,LOS='y'):
    plt.close('all')
    fig,ax = plt.subplots(1,3, sharex=True,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()

    for nf,field in enumerate(['r_TE'+LOS,'r_TB'+LOS,'r_EB'+LOS]):
        for sim in simlist:
            this_sim=simulation.corral[sim]
            this_sim.load()
            qu = field[-3:-1].upper()
            label = r'$r_{%s} \parallel$'%qu
            #axlist[nf  ].text(0.1, -0.8, label)
            label = r'$r_{%s} \perp$'%qu
            axlist[nf].set_title(field)
            #axlist[nf+3].text(0.1, -0.8, label)
            axlist[nf  ].plot(this_sim.avg_spectra['k2d'], this_sim.avg_spectra[field], c=this_sim.color, linestyle=this_sim.linestyle)

            #axlist[nf+3].plot(spectra_dict['y'][sim].lcent, spectra_dict['y'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim])
            #axlist[nf+6].plot(proj.lcent, spectra_dict['z'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim])
    for a in axlist:
        dt.axbonk(a,xscale='log',yscale='log',xlabel=None,ylabel=None)
        #a.set_yscale('symlog',linthresh=0.09)
        a.set_yscale('linear')
        a.set_ylim([-1,1])
        #a.set_ylim([-10,10])
    for a in axlist:
        a.set_xlabel(r'$k/k_max$')
#        a.axhline( 5e-2,c=[0.5]*4)
#        a.axhline(-5e-2,c=[0.5]*4)
#    for a in [axlist[0]]:
#        a.set_ylabel(r'$r_{XY}$')
#        a.set_yticks([-1,-0.1,-0.01,0.01,0.1,1])

    outname = '%s/r_XY.pdf'%dl.plotdir
    fig.savefig(outname)
    print(outname)

#
# reb, ttb, rte PDFs
#

def no():
    plt.close('all')
    fig,ax = plt.subplots(2,3, sharex=True,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()
    bins = np.linspace(-1,1,64)
    linear = 0.1
    b1 = np.linspace(-linear, linear, 8)
    b2 = np.logspace(np.log10(linear), 0, 8)
    bins = np.concatenate([-b2[::-1], b1, b2])
    bins = np.unique(bins)
    for sim in sim_colors.simlist:
        #for a in ax.flatten():
        #    a.clear()
        for nf,field in enumerate(['avg_rte', 'avg_rtb', 'avg_reb']):
            qu = field[-2:].upper()
            label = r'$r_{%s} \parallel$'%qu
            axlist[nf  ].text(-0.5, 100, label)
            label = r'$r_{%s} \perp$'%qu
            axlist[nf+3].text(-0.5, 100, label)
            h1,b1=np.histogram(spectra_dict['x'][sim].spectra[field],bins=bins)
            h2,b2=np.histogram(spectra_dict['y'][sim].spectra[field],bins=bins)
            b1c = 0.5*(b1[:-1]+b1[1:])
            b2c = 0.5*(b2[:-1]+b2[1:])
            b1w = (b1[1:]-b1[:-1])
            b2w = (b2[1:]-b2[:-1])

            axlist[nf  ].plot(b1c, h1, color=sim_colors.color[sim])#, linestyle=sim_colors.linestyle[sim])
            axlist[nf+3].plot(b2c, h2, color=sim_colors.color[sim])#, linestyle=sim_colors.linestyle[sim])
        #plt.savefig("%s/r_XY_pdf_%s.pdf"%(plotdir,sim))
    for a in axlist:
        dt.axbonk(a,yscale='log',xlabel=None,ylabel=None)
        a.set_xscale('symlog',linthresh=0.151515151515151515151515151515)
#        a.set_ylim([-1,1])
        a.set_xlim([-1,1])
    for a in axlist[3:]:
        a.set_xlabel(r'$r_{XY}$')
#   for a in [axlist[0],axlist[3]]:
#       a.set_ylabel(r'$P(r_{XY})$')
#       a.set_yticks([-1,-0.1,-0.01,0.01,0.1,1])

    outname = '%s/h_XY_pdf.pdf'%plotdir
    fig.savefig(outname)
    print(outname)


#
# reb, ttb, rte mean and variance
#

if 0:
    plt.close('all')
    #fig,ax = plt.subplots(2,3, sharex=True,sharey=True)
    fig,ax = plt.subplots(1,3, sharex=True,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()
    for sim in sim_colors.simlist:
        #for a in ax.flatten():
        #    a.clear()
        for nf,field in enumerate(['avg_rte', 'avg_rtb', 'avg_reb']):
            qu = field[-2:].upper()
            label = r'$r_{%s} \parallel$'%qu
            #axlist[nf  ].text(0.5e-3, 0.2, label)
            label = r'$r_{%s} \perp$'%qu
            #axlist[nf+3].text(0.5e-3, 0.2, label)

            mean1 = spectra_dict['x'][sim].spectra[field].mean()
            std1  = spectra_dict['x'][sim].spectra[field].std()
            mean2 = spectra_dict['y'][sim].spectra[field].mean()
            std2  = spectra_dict['y'][sim].spectra[field].std()

            axlist[nf  ].scatter(np.abs(mean1),std1, color=sim_colors.color[sim], marker=sim_colors.marker[sim])
            #axlist[nf+3].scatter(np.abs(mean2),std2, color=sim_colors.color[sim], marker=sim_colors.marker[sim])

    for a in axlist:
        #dt.axbonk(a,xscale='linear',yscale='linear',xlabel=None,ylabel=None)
        dt.axbonk(a,xscale='log',yscale='log',xlabel=None,ylabel=None)
    for na,a in enumerate(axlist):
        val = ['TE','TB','EB'][na]
        a.set_xlabel(r'$\langle r_{%s} \rangle$'%val)
    for a in [axlist[0]]:
        a.set_ylabel(r'$\sigma_{XY}$')

    outname = '%s/r_XY_mean_std.pdf'%plotdir
    fig.tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.savefig(outname)
    print(outname)

