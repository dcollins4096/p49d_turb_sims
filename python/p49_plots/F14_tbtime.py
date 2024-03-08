from GL import *
import tools.davetools as dt
def plot_rtb_time(sim_list):
    plt.close('all')

    nwide=min([3,len(sim_list)])
    ntall = int(np.ceil(len(sim_list)//nwide))
    #for spectra
    f2s, axs = plt.subplots(ntall,nwide, figsize=(12,12))
    axs_list=axs.flatten()
    #for time series
    f2t, axt = plt.subplots(ntall,nwide, figsize=(12,12))
    axt_list=axt.flatten()
    #for histograms
    f2h, axh = plt.subplots(ntall,nwide, figsize=(12,12))
    axh_list=axh.flatten()

    #for the summary plot
    fT, axT = plt.subplots(1,1)


    for ns, sim in enumerate(sim_list):
        print('rtb',sim)
        fig,axes=plt.subplots(1,3)
        ax=axes[0];ax1=axes[1];ax2=axes[2]
        this_sim = simulation.corral[sim]
        xvals = this_sim.avg_spectra['k2d']
        fit_range =this_sim.get_fitrange(xvals)
        mask = (xvals > fit_range[0])*(xvals < fit_range[1])
        rm = dt.rainbow_map(len(this_sim.ann_frames))
        frames = this_sim.quan_time['frames'][this_sim.ann_frame_mask]
        times = this_sim.quan_time['time'][this_sim.ann_frame_mask]
        means = []
        stds = []
        
        all_signal = []
        for nf,frame in enumerate(frames):
            field='r_TBy'
            time = times[nf]
            signal = this_sim.all_spectra[frame][field][mask]
            all_signal += list(signal)
            means.append(signal.mean())
            stds.append(signal.std())

            ax.plot(xvals, this_sim.all_spectra[frame][field], c=rm(nf), linewidth=0.5)
            ax.axvline(fit_range[0],c='k',linewidth=0.2)
            ax.axvline(fit_range[1],c='k',linewidth=0.2)

            ax1.scatter(time,signal[-1], c=[rm(nf)])
            ax1.errorbar(time, signal[-1], yerr=stds[-1], c=rm(nf))

            axs_list[ns].plot(xvals, this_sim.all_spectra[frame][field], c=rm(nf), linewidth=0.5)
            axs_list[ns].axvline(fit_range[0],c='k',linewidth=0.2)
            axs_list[ns].axvline(fit_range[1],c='k',linewidth=0.2)
            axs_list[ns].plot(xvals, this_sim.avg_spectra[field], c='k', linewidth=1)

            axs_list[ns].set(ylim=[-1,1],xscale='log')

            axt_list[ns].scatter(time,signal[-1], c=[rm(nf)])
            axt_list[ns].errorbar(time, signal[-1], yerr=stds[-1], c=rm(nf))

        x = nar(sorted(all_signal))
        y = np.arange(x.size)/x.size
        axh_list[ns].plot(x,y)
        mu = nar(all_signal).mean()
        sig = nar(all_signal).std()

        axT.scatter( this_sim.Ms_mean, mu)
        axT.errorbar( this_sim.Ms_mean, mu, yerr=sig)




        xbins = np.linspace(x.min(),x.max(),128)
        xbins = np.linspace(-10,10,1024)
        dx = xbins[1]-xbins[0]
        g = 1/np.sqrt(2*np.pi*sig**2)*np.exp(-0.5*(xbins-mu)**2/sig**2)
        y2 = np.cumsum(g*dx)
        axh_list[ns].plot(xbins,y2,c='k',linewidth=0.2)
        #axh_list[ns].hist(all_signal,histtype='step')
        axh_list[ns].axhline(0.5,c='k',linewidth=0.2)
        axh_list[ns].axvline(0.0,c='k',linewidth=0.2)
        axh_list[ns].set(xlim=[-0.75,0.75])

        axs_list[ns].set(title='(%.1f,%.1f)'%(this_sim.Ms_nom,this_sim.Ma_nom))
        axt_list[ns].set(title='(%.1f,%.1f)'%(this_sim.Ms_nom,this_sim.Ma_nom), ylim=[-1,1])
        axt_list[ns].axhline(0,c='k',linewidth=0.2)
        axs_list[ns].axhline(0,c='k',linewidth=0.2)




        #ax2.hist(signal,histtype='step')
        ax2.set(xlim=[-1,1])
        ax2.axvline(0)
        ax.set(ylim=[-1,1],xscale='log')
        ax1.set(ylim=[-1,1])
        #fig.savefig('%s/rtb_%s'%(dl.plotdir,sim))
        plt.close(fig)
    f2s.tight_layout()
    f2s.savefig('%s/rtb_allsim.pdf'%(dl.plotdir))
    f2t.tight_layout()
    f2t.savefig('%s/rtb_allsim_time.pdf'%(dl.plotdir))
    f2h.tight_layout()
    f2h.savefig('%s/rtb_allsim_hist.pdf'%(dl.plotdir))
    fT.savefig('%s/rtb_summary.pdf'%(dl.plotdir))



