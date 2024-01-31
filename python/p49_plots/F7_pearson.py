

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
    fig,ax = plt.subplots(1,3,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()

    for nf,field in enumerate(['r_TE'+LOS,'r_TB'+LOS,'r_EB'+LOS]):
        collector=[]
        for sim in simlist:
            this_sim=simulation.corral[sim]
            this_sim.load()
            qu = field[-3:-1].upper()
            label = r'$r_{%s}\ \hat{%s}$'%(qu, LOS)
            axlist[nf].set_title(label)
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
            #axlist[nf].scatter(np.abs(mean),std, c=[this_sim.color],marker=this_sim.marker,edgecolor=ec, s=this_sim.marker_size*20)
            axlist[nf].scatter(mean,std, c=[this_sim.color],marker=this_sim.marker,edgecolor=ec, s=this_sim.marker_size*20)
            axlist[nf].set(xlabel=r'$\mu$', ylabel=None)
            if nf==0:
                axlist[nf].axvline(0.35,c=[0.5]*3)

            collector.append(mean)

        if 1:
            j = nar(collector)
            j.sort()
            print(field,'mode',j[j.size//2])
            y = np.arange(j.size)/j.size
            twin=axlist[nf].twinx()
            twin.plot(j,y,c='r')
            #twin.axhline(0.5)
            twin.set(ylim=[0,1], yticks=[])
            
            #axlist[nf+3].plot(spectra_dict['y'][sim].lcent, spectra_dict['y'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim])
            #axlist[nf+6].plot(proj.lcent, spectra_dict['z'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim])
    twin.set(yticks=np.linspace(0,1,11), ylabel=r'$Cuml(\mu)$')
    for naa,a in enumerate(axlist):
        #dt.axbonk(a,xscale='log',yscale='log',xlabel=None,ylabel=None)
        #a.set_yscale('symlog',linthresh=0.09)
        a.set_yscale('linear')
        #a.set_ylim([-1,1])
        #a.set_ylim([-10,10])
    axlist[0].set(ylabel=r'$\sigma_{XY}$')
    axlist[0].set(xscale='linear')
    #kludge
    #axlist[1].set(xscale='linear', xlim=[0,0.05])
    #axlist[2].set(xscale='linear', xlim=[0,0.05])
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
            label = r'$r_{%s}\ \hat{%s}$'%(qu, LOS)
            axlist[nf].set_title(label)
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
    fig,axes = plt.subplots(1,3,figsize=(8,4))
    axlist=axes.flatten()

    for nf,field in enumerate(['r_TE'+LOS,'r_TB'+LOS,'r_EB'+LOS]):
        for sim in simlist:
            this_sim=simulation.corral[sim]
            this_sim.load()
            qu = field[-3:-1].upper()
            label = r'$r_{%s}\ \hat{%s}$'%(qu, LOS)
            label = r'$r_{%s}$'%(qu)
            axlist[nf].set_ylabel(label)
            #axlist[nf+3].text(0.1, -0.8, label)
            xvals=this_sim.avg_spectra['k2d']
            axlist[nf  ].plot(xvals, this_sim.avg_spectra[field], c=this_sim.color, linestyle=this_sim.linestyle)

            fitrange = this_sim.get_fitrange(xvals)

    for a in axlist:
        a.set(xscale='log',yscale='log')
        #a.set_yscale('symlog',linthresh=0.09)
        a.set_yscale('linear')
        a.set_ylim([-.25,.25])
        a.axhline(0,c=[0.5]*4)
        a.axvline(fitrange[0], c=[0.5]*4)
        a.axvline(fitrange[1], c=[0.5]*4)
        #a.set_ylim([-10,10])
    axlist[0].set_ylim([-1,1])
    for a in axlist:
        a.set_xlabel(r'$k$')
#        a.axhline( 5e-2,c=[0.5]*4)
#        a.axhline(-5e-2,c=[0.5]*4)
#    for a in [axlist[0]]:
#        a.set_ylabel(r'$r_{XY}$')
#        a.set_yticks([-1,-0.1,-0.01,0.01,0.1,1])

    outname = '%s/r_XY.pdf'%dl.plotdir
    fig.tight_layout()
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

def plot_machmean(simlist,LOS='y'):
    plt.close('all')
    fig,ax = plt.subplots(1,3, figsize=(8,3))
    #fig.subplots_adjust(wspace=0, hspace=0)

    axlist=ax.flatten()

    ext=dt.extents()
    from collections import defaultdict
    rte_collector=defaultdict(list)
    for nf,field in enumerate(['r_TE'+LOS,'r_TB'+LOS,'r_EB'+LOS]):
        collector=[]
        for sim in simlist:
            this_sim=simulation.corral[sim]
            this_sim.load()
            qu = field[-3:-1].upper()
            #label = r'$r_{%s}\ \hat{%s}$'%(qu, LOS)
            label = r'$\langle r_{%s}\rangle$'%(qu)
            #label = r'$r_{%s}$'%(qu)
            xvals = this_sim.avg_spectra['k2d']
            fit_range =this_sim.get_fitrange(xvals)
            mask = (xvals > fit_range[0])*(xvals < fit_range[1])
            signal = this_sim.avg_spectra[field][mask]
            mean = signal.mean()
            std  = signal.std()
            if mean < 0:
                ec=None
            else:
                ec=None
            #axlist[nf].scatter(np.abs(mean),std, c=[this_sim.color],marker=this_sim.marker,edgecolor=ec, s=this_sim.marker_size*20)
            axlist[nf].scatter(this_sim.Ms_mean,mean, c=[this_sim.color],marker=this_sim.marker,edgecolor=ec, s=this_sim.marker_size*20)
            axlist[nf].errorbar(this_sim.Ms_mean,mean, yerr=std, c=this_sim.color)
            axlist[nf].set(xlabel=sim_colors.mach_label, ylabel=label)
            rte_collector[field].append(mean)

            if nf>0:
                ext(mean)
            collector.append(mean)

    for f in rte_collector:
        print(f)
        rte=nar(rte_collector[f])
        print(rte.mean())
        print(rte.std())
    for nf in [0,1,2]:
        fid_line=[0.355,0.05,None][nf]
        if fid_line is not None:
            axlist[nf].axhline(fid_line,c=[0.5]*4)

    for naa,a in enumerate(axlist):
        #dt.axbonk(a,xscale='log',yscale='log',xlabel=None,ylabel=None)
        #a.set_yscale('symlog',linthresh=0.09)
        a.set_yscale('linear')
        #a.set_ylim([-1,1])
        #a.set_ylim([-10,10])
        if naa>0:
            a.set(ylim=ext.minmax)
    axlist[0].set(xscale='linear')
    #kludge
    #axlist[1].set(xscale='linear', xlim=[0,0.05])
    #axlist[2].set(xscale='linear', xlim=[0,0.05])
#    for a in axlist:
#        a.set_xlabel(r'$r_{XY}$')
#        a.axhline( 5e-2,c=[0.5]*4)
#        a.axhline(-5e-2,c=[0.5]*4)
#    for a in [axlist[0]]:
#        a.set_ylabel(r'$r_{XY}$')
#        a.set_yticks([-1,-0.1,-0.01,0.01,0.1,1])

    outname = '%s/mach_mean_rXY_%s.pdf'%(dl.plotdir,LOS)
    fig.tight_layout()
    fig.savefig(outname)
    print(outname)

def plot_tb_eb(simlist,LOS='y', ax=None):

    save_fig=False
    if ax is None:
        fig,ax=plt.subplots(1,1)
        save_fig=True
    for sim in simlist:
        this_sim=simulation.corral[sim]
        this_sim.load()
        xvals = this_sim.avg_spectra['k2d']
        fit_range =this_sim.get_fitrange(xvals)
        mask = (xvals > fit_range[0])*(xvals < fit_range[1])
        TB = this_sim.avg_spectra['r_TB'+LOS][mask].mean()
        EB = this_sim.avg_spectra['r_EB'+LOS][mask].mean()
        ax.scatter(TB, EB, c=[this_sim.color],marker=this_sim.marker, s=this_sim.marker_size*20)
    ax.set(xlabel=r'$r^{TB}$',ylabel = r'$r^{EB}$')
    colorbar=plt.colorbar(sim_colors.cbar,ax=ax)
    colorbar.set_label(sim_colors.mach_label)
    if save_fig:
        fig.savefig('%s/TB_EB.pdf'%dl.plotdir)

