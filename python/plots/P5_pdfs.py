
from GL import *
import simulation
root4pi = np.sqrt(4*np.pi)

def plot_sigma(simlist):

    fig,axes=plt.subplots(1,2, figsize=(8,4))
    ax=axes
    ax0=axes[0];ax1=axes[1]#;ax2=axes[2]
    
    fields=['magnetic_field_%s'%s for s in 'xyz']

    for sim in simlist:
        this_sim=simulation.corral[sim]
        this_sim.load()
        #this_sim.read_pdfs(fields, pdf_prefix=pdf_prefix)
        mask = this_sim.ann_frame_mask
        sbx = this_sim.quan_time['bx_std'][mask]
        sby = this_sim.quan_time['by_std'][mask]
        sbz = this_sim.quan_time['bz_std'][mask]
        vrms= this_sim.quan_time['vrms'][mask]
        kwargs={'color':this_sim.color,'marker':this_sim.marker}
        x=nar([0.6,1.9])
        y1=x
        y2=x-0.1
        y3=x+0.1
        ax0.plot(x,y1,c='k')
        ax0.plot(x,y2,c='k')
        ax0.plot(x,y3,c='k')
        ax0.scatter(sby/sbx,sbz/sbx, **kwargs)
        ax1.scatter(vrms, sby/sbz, **kwargs)
    ax0.set(xlabel=r'$\sigma_{by}/\sigma_{bx}$',ylabel=r'$\sigma_{bz}/\sigma_{bx}$')
    ax1.set(xlabel=r'$v_{rms}$',ylabel=r'$\sigma_{by}/\sigma_{bz}$')
    fig.savefig('%s/magnetic_sigmas'%(dl.plotdir))


    


def plot_pdfs(simlist,fields, name="STUFF", pdf_prefix="pdf_scaled", all_or_ann_frames='all', norm_axis=False, overgauss=False,
              logy=False,plot_all=True, all_on_one=False):

    total_columns = min([3,len(simlist)])
    total_rows = len(simlist)//total_columns
    wrong_bmag_norm=True


    fig,axes=plt.subplots(total_rows, total_columns, figsize=(20,20))
    if hasattr(axes,'size'):
        axlist=axes.flatten()
    else:
        axlist=[axes]
    ext_x=dt.extents()
    for nfield, field in enumerate(fields):
        if not all_on_one:
            for ax in axlist:
                ax.clear()

        print(field, fields)
        for nsim, sim in enumerate(simlist):
            this_sim = simulation.corral[sim]
            this_sim.load()

            this_sim.read_pdfs(fields, pdf_prefix=pdf_prefix)
            scaled=False
            if pdf_prefix in ['pdf_scaled']:
                scaled=True

            thax=axlist[nsim]

            def puller(arr):
                return arr[this_sim.frame_mask].mean()
            def puller2(arr):
                return np.sqrt( (arr[this_sim.frame_mask]**2).mean())
            measured_mean = this_sim.quan_mean['bx_avg']
            sbx,sby,sbz=[this_sim.quan_mean['b%s_std'%s] for s in 'xyz']
            measured_sigma=np.mean([sbx,sby,sbz])

            local_average = 0
            npdf=0
            if all_or_ann_frames=='all':
                frames=this_sim.pdfs[field].keys()
            else:
                frames=this_sim.ann_frames
                #frames=frames[-1:]


            for frame in frames:
                hist = this_sim.pdfs[field][frame]['hist']
                cbins= this_sim.pdfs[field][frame]['cbins']
                local_mean = (cbins*hist).sum()/hist.sum()
                if hist.sum()==0:
                    print("Bad histogram", sim, frame)
                #if np.abs(local_mean) > 1e-3:
                #    print("Bad Frame",sim,frame,local_mean)
                #print("%10d %0.2f"%(frame,local_mean))
                local_average = hist + local_average
                npdf+=1
                the_x=cbins
                if norm_axis:
                    the_x=cbins/measured_sigma
                if wrong_bmag_norm and field == 'magnetic_field_strength':
                    WRONG=np.sqrt(3)
                    the_x=cbins#*WRONG

                the_y=hist/hist.max()
                ok = slice(None)
                if overgauss:
                    my_mean=0
                    my_sigma=1
                    GGG = np.exp(-(cbins-my_mean)**2/(2*my_sigma**2))
                    the_y/=GGG
                    if logy:
                        thax.set(yscale='log', xlim=[-5,5], ylim=[1e-3,1e3])
                    ok = the_y>0
                if plot_all:
                    thax.plot(the_x[ok],the_y[ok], c=[0.5]*3, linewidth=0.5)
                ext_x(the_x)

            thax.text(2,0.8,this_sim.name)
            #avg_pdf=this_sim.avg_pdf[field]
            avg_pdf=local_average/npdf
            the_y = avg_pdf/avg_pdf.max()
            if overgauss:
                the_y/=GGG
                thax.axhline(2,c='k')
                thax.axhline(0.5,c='k')
            thax.plot(the_x,the_y, c='r', linewidth=1)
            pdf_mean = (cbins*avg_pdf).sum()/avg_pdf.sum()
            pdf_sigma = np.sqrt(( (cbins-pdf_mean)**2*avg_pdf).sum()/avg_pdf.sum())
            #print("PDF mean vs Avg mean", mean/ bx)
            if 0:
                #come back to this.
                print("PDF sigma vs avg", measured_sigma/pdf_sigma)

            if field in ['magnetic_field_x','magnetic_field_y','magnetic_field_z'] and not overgauss:
                #my_sigma=pdf_sigma
                #my_mean =pdf_mean
                #my_sigma=measured_sigma
                #my_mean=measured_mean
                if scaled:
                    my_mean=0
                    my_sigma=1
                GGG = np.exp(-(cbins-my_mean)**2/(2*my_sigma**2))
                thax.plot(the_x,GGG/GGG.max(),c='g', linewidth=1)
                #thax.set(yscale='log', ylim=[1e-5,1.01])
            if field in ['magnetic_field_strength'] and False:
                #Dist = cbins*np.exp(+mu_bx/sigma_bx1*cbins-cbins**2/(2*sigma_bx1))
                #The sinh blows up, so keep it real.
                arg = measured_mean*cbins/(measured_sigma)
                ok = arg<500
                jbins=cbins[ok]
                arg = measured_mean*jbins/(measured_sigma)

                Dist = jbins*np.exp( -jbins**2/(2*measured_sigma**2))*np.sinh(arg)
                this_sigma = measured_sigma
                mean = measured_mean
                if wrong_bmag_norm:
                    WRONG=np.sqrt(3)
                    the_x_in = cbins/WRONG
                    the_x = cbins/WRONG
                else:
                    WRONG=1
                    the_x_in = cbins
                #Thing = cbins*np.exp( -cbins**2/(2*this_sigma))*np.sinh(mean*cbins/(this_sigma))
                if scaled:
                    Thing = cbins*np.exp( -the_x_in**2/(2))*np.sinh(mean*cbins/(this_sigma*WRONG))
                    Thing2 = cbins*np.exp( -WRONG**2*the_x_in**2/(2)*this_sigma)*np.sinh(pdf_mean*cbins*WRONG)
                else:
                    Thing = cbins*np.exp( -cbins**2/(2*this_sigma**2))*np.sinh(mean*cbins/(this_sigma**2))
                    Thing2 = cbins*np.exp( -cbins**2/(2*this_sigma))*np.sinh(pdf_mean*cbins/(this_sigma))
                    Thing3 = cbins*np.exp( -the_x**2/2*this_sigma)*np.sinh(pdf_mean*the_x)

                Dist=Thing
                thax.plot(the_x,Dist/Dist.max(),c='g',linewidth=0.5)
                thax.plot(the_x,Thing2/Thing2.max(),c='purple',linewidth=0.5)
                #thax.plot(the_x,Thing3/Thing3.max(),c='k',linewidth=1)

            thax.set(xlabel=field)
        outname = "%s/pdf_%s_%s"%(dl.plotdir,name,field)

        #for ax in axes.flatten():
        #    #ax.set_xlim(ext_x.minmax)
        #    ax.set(xlim=[0,12])

        
        fig.savefig(outname)
        print(outname)
        plt.close(fig)

def plot_pdfs_fits(simlist, name="STUFF",norm_axis=False):
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
                #don't do it this way, use quan.
                #UNITS=root4pi
                #bx,by,bz =    [UNITS*last(this_sim.quan_time['b%s_avg'%s]) for s in 'xyz']
                #sbx,sby,sbz = [UNITS*last(this_sim.quan_time['b%s_std'%s]) for s in 'xyz']
                #vx,vy,vz =    [last(this_sim.quan_time['v%s_std'%s]) for s in 'xyz']
                #v3d = np.sqrt(vx*vx+vy*vy+vz*vz)
                #b3d = np.sqrt(sbx**2+sby**2+sbz**2)
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
                    measured_mean = this_sim.quan_mean['bx_avg']
                    sbx,sby,sbz=[this_sim.quan_mean['b%s_std'%s] for s in 'xyz']
                    measured_sigma=np.mean([sbx,sby,sbz])

                    this_sigma=sum(sigmas)/3
                    #print(this_sigma/measured_sigma)
                    #print("mean, meas %0.2f %0.2f Bnom %0.2f"%(mean,measured_mean,this_sim.B_nom))
                    #expect = 0.5*(measured_mean+np.sqrt(measured_mean**2-4*measured_sigma**2))
                    #print("mean, pdf/exp %0.2f meas/nom %0.2f"%(mean/expect,measured_mean/this_sim.B_nom))
                    Max = cbins**2*np.exp( -(cbins-mean)**2/(2*this_sigma**2))

                    #the one the fits ok
                    Thing = cbins*np.exp( -cbins**2/(2*this_sigma))*np.sinh(mean*cbins/(this_sigma))
                    #the "right" one.
                    BetterThing = cbins*np.exp( -cbins**2/(2*measured_sigma**2))*np.sinh(measured_mean*cbins/(measured_sigma**2))
                    #BetterThing = cbins*np.exp( -cbins**2/(2*measured_sigma))*np.sinh(mean*cbins/(measured_sigma))

                    if norm_axis:
                        nrm=this_sigma
                    else:
                        nrm=1
                    #nrm=1
                    #thax.axvline(bx)
                    the_x=cbins/nrm
                    thax.plot(the_x,mhist, c=this_sim.color,linestyle=this_sim.linestyle)
                    thax.text(2,0.8,this_sim.name)
                    ext(the_x)
                    thax.plot(the_x,BetterThing/BetterThing.max(),c='k')

                    ax2.scatter(mean, this_sigma, **kwargs)
                    ax2b.scatter(this_sim.Ms_nom,np.sqrt(4*np.pi)*this_sim.Ms_nom/mean/this_sim.Ma_nom, **kwargs)
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
