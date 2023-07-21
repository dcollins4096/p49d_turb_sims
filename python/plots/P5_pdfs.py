
from GL import *
import simulation
root4pi = np.sqrt(4*np.pi)

def plot_pdfs(simlist,fields, name="STUFF", pdf_prefix="pdf_scaled"):

    total_columns = min([3,len(simlist)])
    total_rows = len(simlist)//total_columns


    for nfield, field in enumerate(fields):
        fig,axes=plt.subplots(total_rows, total_columns, figsize=(20,20))
        if hasattr(axes,'size'):
            axlist=axes.flatten()
        else:
            axlist=[axes]
        for nsim, sim in enumerate(simlist):
            this_sim = simulation.corral[sim]
            this_sim.load()

            this_sim.read_pdfs(fields, pdf_prefix=pdf_prefix)

            thax=axlist[nsim]

            def puller(arr):
                return arr[this_sim.frame_mask].mean()
            def puller2(arr):
                return np.sqrt( (arr[this_sim.frame_mask]**2).mean())
            UNITS=root4pi
            bx,by,bz =    [UNITS*puller(this_sim.quan_time['b%s_avg'%s]) for s in 'xyz']
            sbx,sby,sbz = [UNITS*puller2(this_sim.quan_time['b%s_std'%s]) for s in 'xyz']
            vx,vy,vz =    [puller(this_sim.quan_time['v%s_std'%s]) for s in 'xyz']
            v3d = np.sqrt(vx*vx+vy*vy+vz*vz)
            b3d = np.sqrt(sbx**2+sby**2+sbz**2)
            mean = bx
            this_sigma=np.mean([sbx,sby,sbz])

            mean_spectra = 0
            for frame in this_sim.pdfs[field]:
                hist = this_sim.pdfs[field][frame]['hist']
                cbins= this_sim.pdfs[field][frame]['cbins']

                thax.plot(cbins,hist/hist.max(), c=[0.5]*3, linewidth=0.5)
            thax.text(2,0.8,this_sim.name)
            avg_pdf=this_sim.avg_pdf[field]
            thax.plot(cbins,avg_pdf/avg_pdf.max(), c='r', linewidth=0.5)
            mean = (cbins*avg_pdf).sum()/avg_pdf.sum()
            print(mean/ bx)
            this_sigma = np.sqrt(( (cbins-mean)**2*avg_pdf).sum()/avg_pdf.sum())

            if field in ['magnetic_field_x','magnetic_field_y','magnetic_field_z']:
                GGG = np.exp(-cbins**2/2)
                thax.plot(cbins,GGG/GGG.max(),c='g', linewidth=0.5)
            if field in ['magnetic_field_strength']:
                #Dist = cbins*np.exp(+mu_bx/sigma_bx1*cbins-cbins**2/(2*sigma_bx1))
                arg = mean*cbins/(this_sigma)
                ok = arg<500
                jbins=cbins[ok]
                arg = mean*jbins/(this_sigma)

                Dist = jbins*np.exp( -jbins**2/(2*this_sigma))*np.sinh(arg)
                thax.plot(cbins[ok],Dist/Dist.max(),c='g',linewidth=0.5)

            thax.set(xlabel=field)
        outname = "%s/pdf_%s_%s"%(dl.plotdir,name,field)
        fig.savefig(outname)
        print(outname)
        plt.close(fig)

def plot_pdfs_fits(simlist, name="STUFF"):
    fields =['magnetic_field_x','magnetic_field_y','magnetic_field_z','magnetic_field_strength' ]

    if 0:
        total_columns = len(fields)
        total_rows = len(simlist)
    total_columns=3
    total_rows=len(simlist)//total_columns

    fig,axes=plt.subplots(total_rows, total_columns, figsize=(20,20))
    axlist=axes.flatten()
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
        for nfield, field in enumerate(fields):
            #if total_rows > 1:
            #    thax = axes[nsim][nfield]
            #else:
            #    thax = axes[nfield]
            thax = axes[nsim//total_columns][ nsim%total_columns]

            kwargs={'c':this_sim.color,'marker':this_sim.marker}
            for frame in this_sim.pdfs[field]:
                hist = this_sim.pdfs[field][frame]['hist']
                cbins= this_sim.pdfs[field][frame]['cbins']
                mean = (cbins*hist).sum()/hist.sum()
                sigma = np.sqrt(( (cbins-mean)**2*hist).sum()/hist.sum())
                sigmas.append(sigma)
                means = max([mean,means])
                mhist = hist/hist.max()
                #thax.plot((cbins-mean)/sigma,mhist, c=this_sim.color,linestyle=this_sim.linestyle)

                avg_index = np.where(this_sim.quan_time['frames'] == frame)[0][0]
                def puller(arr):
                    return arr[this_sim.frame_mask].mean()
                def puller2(arr):
                    return np.sqrt( (arr[this_sim.frame_mask]**2).mean())
                def last(arr):
                    return arr[-1]
                UNITS=root4pi
                bx,by,bz =    [UNITS*last(this_sim.quan_time['b%s_avg'%s]) for s in 'xyz']
                sbx,sby,sbz = [UNITS*last(this_sim.quan_time['b%s_std'%s]) for s in 'xyz']
                vx,vy,vz =    [last(this_sim.quan_time['v%s_std'%s]) for s in 'xyz']
                v3d = np.sqrt(vx*vx+vy*vy+vz*vz)
                b3d = np.sqrt(sbx**2+sby**2+sbz**2)
                if 0:
                    #
                    # what is going on with means?
                    # future me broke the b3d code, need to take the root4pi back out
                    #
                    if field[-1] == 'x' and False:
                        Gaus = np.exp( -(cbins-mean)**2/(2*sigma**2))
                        thax.plot( cbins, Gaus/Gaus.max())
                        thax.axvline(bx*root4pi, c='r')
                        thax.axvline(mean,c='magenta')
                        ext(cbins)

                    if field[-1]=='x' and True:
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

            if 1:
                if field[-1] in 'xyz':
                    #don't make the xyz pdfs at this time.
                    pass

                    #Gaus = np.exp( -(cbins-mean)**2/(2*sigma**2))
                    #thax.plot( cbins, Gaus/Gaus.max())
                else:
                    this_sigma=sum(sigmas)/3
                    #print(this_sigma/((sbx+sby+sbz)/3))
                    Max = cbins**2*np.exp( -(cbins-mean)**2/(2*this_sigma**2))

                    Thing = cbins*np.exp( -cbins**2/(2*this_sigma))*np.sinh(mean*cbins/(this_sigma))
                    #Thing_w1 = np.exp( (2*mean*cbins-cbins**2)/(2*this_sigma))
                    #Thing_w1 = np.exp( (-(cbins-mean)**2)/(2*this_sigma))
                    #Thing_w2 = cbins*np.exp( -cbins**2/(2*this_sigma))*np.exp(mean*cbins/(this_sigma))
                    #Thing_w3 = cbins**2*np.exp( -cbins**2/(2*this_sigma))#*np.exp(mean*cbins/(this_sigma))

                    #thax.plot(cbins,Max/Max.max(), color='magenta')
                    #thax.plot(cbins,Max2/Max2.max(), color='magenta', linestyle=":")
                    nrm=this_sigma
                    nrm=this_sigma
                    #thax.axvline(bx)
                    the_x=cbins/nrm
                    thax.plot(the_x,mhist, c=this_sim.color,linestyle=this_sim.linestyle)
                    thax.text(2,0.8,this_sim.name)
                    ext(the_x)

                    ax2.scatter(mean, this_sigma, **kwargs)
                    ax2b.scatter(this_sim.ms,np.sqrt(4*np.pi)*this_sim.ms/mean/this_sim.ma, **kwargs)
                    #print(sim,mean)

                    thax.plot(the_x,Thing/Thing.max(),color='cyan')
                    #thax.plot(the_x,Thing_w1/Thing_w1.max(),color='cyan')
                    #thax.plot(the_x,Thing_w2/Thing_w2.max(),color='cyan')
                    #thax.plot(the_x,Thing_w3/Thing_w3.max(),color='cyan')

                    #print(bx,by,bz)
                    #print(sim,bx/mean*root4pi)

                #thax.plot(cbins,hist)

            thax.set(xlabel=field)

    fig2.tight_layout()
    fig2.savefig('plots_to_sort/mean_sigma.png')
    for ax in axlist:
        ax.set_xlim(ext.minmax)
    outname = "%s/pdf_%s"%(dl.plotdir,name)
    fig.tight_layout()
    fig.savefig(outname)
    print(outname)

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
