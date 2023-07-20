
from GL import *
import simulation
root4pi = np.sqrt(4*np.pi)


def mean_var(simlist, name="STUFF"):
    fields =['magnetic_field_x','magnetic_field_y','magnetic_field_z','magnetic_field_strength' ]

    fig2,axes2=plt.subplots(2,2)
    ax2=axes2[0][0]
    ax2b=axes2[0][1]
    ax2c=axes2[1][0]
    ax2d=axes2[1][1]

    ext=dt.extents()
    for nsim, sim in enumerate(simlist):
        #if sim[-1]!= 'f':
        #    continue
        print(sim)
        this_sim = simulation.corral[sim]
        this_sim.load()

        this_sim.read_pdfs(fields, pdf_prefix='pdf')

        sigmas=[]
        means = 0

        kwargs={'c':this_sim.color,'marker':this_sim.marker}

        def puller(arr):
            return arr[this_sim.frame_mask].mean()
        def puller2(arr):
            return np.sqrt( (arr[this_sim.frame_mask]**2).mean())
        bx,by,bz =    [puller(this_sim.quan_time['b%s_avg'%s]) for s in 'xyz']
        sbx,sby,sbz = [puller2(this_sim.quan_time['b%s_std'%s]) for s in 'xyz']
        vx,vy,vz =    [puller(this_sim.quan_time['v%s_std'%s]) for s in 'xyz']
        v3d = np.sqrt(vx*vx+vy*vy+vz*vz)
        b3d = np.sqrt(sbx**2+sby**2+sbz**2)
        #Plot sigma_b vs sigma_v and mean_b in several ways.
        #ax2.scatter(mean, 1-mean/(bx*root4pi))
        #ax2b.scatter(this_sim.ms,this_sim.ms/v3d,**kwargs)
        #ax2.plot(np.arange(0.5,6),np.arange(0.5,6))
        ax2d.scatter(v3d, bx*b3d/v3d**2, **kwargs)
        ax2d.set(ylim=[0,1])
        ax2d.set(xlabel='vrms',ylabel=r'$\mu \sigma_{bx}/vrms^2$')
        #ax2.plot(np.arange(0.5,7),np.arange(0.5,7),c='k')

        ax2.scatter(v3d, b3d, **kwargs)
        ax2.set(xlabel='vrms',ylabel='brms')
        ax2b.scatter(v3d, b3d/bx,**kwargs)
        ax2b.set(xlabel='vrms',ylabel='sigma/mu')
        ax2c.scatter(v3d, (b3d*bx)/v3d,**kwargs)
        ax2c.set(xlabel='vrms',ylabel=r'$\mu \sigma/vrms$')


    fig2.tight_layout()
    fig2.savefig('plots_to_sort/mean_sigma.png')

def stack_pdfs_raw(simlist,fields, name="STUFF"):
    #
    # Mostly code for making sigma_b and mu_b array.  
    # Don't break it.
    # Final frame PDFs are in pdf_magnetic_field_xyz.h5
    #
    #

    total_columns = len(fields)
    total_rows = 1

    fig,axes=plt.subplots(total_rows, total_columns, figsize=(12,8))
    fig2,axes2=plt.subplots(2,3)

    fid_bins = np.linspace(-10,10,256)
    fid_gauss = np.exp(-fid_bins**2/2)

    fid_sigma={}
    fid_means = {}
    for nsim, sim in enumerate(simlist):
        fid_sigma[sim]=0
        this_sim = simulation.corral[sim]
        this_sim.load()

        this_sim.read_pdfs(fields, pdf_prefix='pdf')

        sigmas={}
        for nfield, field in enumerate(fields):
            thax = axes[nfield]
            thax.plot(fid_bins,fid_gauss,c='k')

            for frame in this_sim.pdfs[field]:
                hist = this_sim.pdfs[field][frame]['hist']
                cbins= this_sim.pdfs[field][frame]['cbins']
                mean = (cbins*hist).sum()/hist.sum()
                sigma = np.sqrt(( (cbins-mean)**2*hist).sum()/hist.sum())
                fid_sigma[sim]+=sigma
                if field == 'magnetic_field_x':
                    fid_means[sim]=mean
                sigmas[field]=sigma
                mhist = hist/hist.max()
                thax.plot((cbins-mean)/sigma,mhist, c=this_sim.color,linestyle=this_sim.linestyle)
                kwargs = {'color':this_sim.color, 'marker':this_sim.marker}
                axes2[0][0].scatter( this_sim.ms, mean,**kwargs)
                axes2[1][0].scatter( this_sim.ma, mean,**kwargs)
                axes2[0][1].scatter( this_sim.ms, sigma,**kwargs)
                axes2[1][1].scatter( this_sim.ma, sigma,**kwargs)
                thax.set(xlabel=field, yscale='log')
                #thax.plot((cbins),hist)
        axes2[0][2].scatter( sigmas['magnetic_field_y']/sigmas['magnetic_field_x'], sigmas['magnetic_field_z']/sigmas['magnetic_field_x'], **kwargs)
        axes2[0][2].plot([0.5,1.5],[0.5,1.5],c='k')

        axes2[0][0].set(xlabel='Ms',ylabel=r'$\mu$')
        axes2[1][0].set(xlabel='Ma',ylabel=r'$\mu$')
        axes2[0][1].set(xlabel='Ms',ylabel=r'$\sigma$')
        axes2[1][1].set(xlabel='Ma',ylabel=r'$\sigma$')
    outname = "%s/pdf_%s.pdf"%(dl.plotdir,name)
    fig.savefig(outname)
    print(outname)
    print('wtf')
    outname = "%s/sigmas_%s.pdf"%(dl.plotdir,name)
    fig2.savefig(outname)
    print(outname)
    if 0:
        for sim in simlist:
            print("sigma_b['%s']=%8.4f"%(sim, fid_sigma[sim]/3))
        for sim in simlist:
            print("mu_bx['%s']=%8.4f"%(sim,fid_means[sim]))
    if 1:
        shooft=[]
        for sim in simlist:
            shooft.append( fid_means[sim]/np.sqrt(2*fid_sigma[sim]))
        print(sorted(shooft))

def pdf_totals(simlist, name="STUFF"):
    fields=['magnetic_field_strength']

    total_columns = 3
    total_rows = len(simlist)//3

    fig,axes=plt.subplots(total_rows, total_columns)

    for nsim, sim in enumerate(simlist):
        this_sim = simulation.corral[sim]
        this_sim.load()
        this_sim.read_pdfs(fields, pdf_prefix='pdf')
        nrow = nsim//total_columns
        ncol = nsim%total_columns
        thax=axes[nrow][ncol]
        field = 'magnetic_field_strength'

        for frame in this_sim.pdfs[field]:
            hist = this_sim.pdfs[field][frame]['hist']
            cbins= this_sim.pdfs[field][frame]['cbins']
            mean = (cbins*hist).sum()/hist.sum()
            sigma = np.sqrt(( (cbins-mean)**2*hist).sum()/hist.sum())
            maxwell = cbins**2*np.exp(-cbins**2/(2*sigma**2))

            thax.plot(cbins,hist/hist.max(),color='magenta')
            think = cbins*np.exp(-cbins**2/(2*sigma**2))*np.exp(mean*cbins/(2*sigma**2))
            thax.plot(cbins,maxwell/maxwell.max(),color='b')
            thax.plot(cbins,think/think.max(),c='g')

        thax.set(xlabel=field)
    outname = "%s/pdf_%s"%(dl.plotdir,name)
    fig.savefig(outname)
    print(outname)
