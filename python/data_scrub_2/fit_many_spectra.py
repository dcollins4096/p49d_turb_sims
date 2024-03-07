from GL import *
from queb3 import powerlaw_fit as plfit

import simulation

def fitter(simlist, lefts=None,rights=None):

    for nsim,sim in enumerate(simlist):
        this_sim=simulation.corral[sim]
        print("Fit!",sim)
        this_sim.load()
        fit_many_spectra(this_sim, lefts, rights)

def fit_many_spectra(self, lefts, rights):

    self.slope_lefts=nar(lefts)
    self.slope_rights=nar(rights)
    self.slope_tensor={}
    for field in self.products_positive:
        self.slope_tensor[field] = np.zeros([len(self.ann_frames), len(lefts), len(rights)])

        for nF, frame in enumerate(self.ann_frames):
            k3d = self.all_spectra[frame]['k3d']
            k2d = self.all_spectra[frame]['k2d']
            if field in self.products_3d:
                xvals = k3d
            else:
                xvals = k2d
            for nL, L in enumerate(self.slope_lefts):
                for nR, R in enumerate(self.slope_rights):

                    spec = self.all_spectra[frame][field]
                    fitrange = [xvals[L],xvals[R]]
                    slope,amp,res=plfit(xvals,spec,fitrange)
                    self.slope_tensor[field][nF,nL,nR]=slope

    fname = "%s/slope_tensor.h5"%self.product_location
    h5ptr=h5py.File(fname,'w')
    try:
        h5ptr['lefts']=self.slope_lefts
        h5ptr['rights']=self.slope_rights
        h5ptr['frames']=self.ann_frames
        for field in self.slope_tensor:
            h5ptr[field]=self.slope_tensor[field]
        print("Write %s"%fname)
    except:
        raise
    finally:
        h5ptr.close()

def plot_spectra(simlist, what='hist', field=None):

    total_columns = min([3,len(simlist)])
    total_rows = len(simlist)//total_columns
    fig,axes=plt.subplots(total_rows, total_columns, figsize=(20,20))
    if hasattr(axes,'size'):
        axlist=axes.flatten()
    else:
        axlist=[axes]

    for nsim,sim in enumerate(simlist):
        this_sim=simulation.corral[sim]
        this_sim.load_spectra_tensor()
        ax = axlist[nsim]
        if what=='hist':
            plotname='hist'
            F = this_sim.slope_tensor[field].flatten()
            F_mean=F.mean()
            F_std=F.std()

            hist,bins,patches=ax.hist((F-F_mean)/F_std, histtype='step')
            ax.axvline(0,c=[0.5]*4)
        if what=='frame':
            plotname='frame'
            F = this_sim.slope_tensor[field]
            Fbar = F.mean(axis=1).mean(axis=1)
            Fbar.shape=Fbar.size,1,1

            std = np.sqrt( ((F-Fbar)**2).mean())
            X = this_sim.ann_frames
            ax.errorbar( X, Fbar, yerr=std)
        if what=='lefts':
            plotname='lefts'
            F = this_sim.slope_tensor[field]
            Fbar = F.mean(axis=0).mean(axis=1)
            Fbar.shape=1,Fbar.size,1
            std = np.sqrt( ((F-Fbar)**2).mean())
            X = this_sim.slope_lefts
            Fbar.shape=Fbar.size
            ax.errorbar( X, Fbar, yerr=std)

        if what=='rights':
            plotname='rights'
            F = this_sim.slope_tensor[field]
            Fbar = F.mean(axis=0).mean(axis=0)
            Fbar.shape=1,1,Fbar.size
            std = np.sqrt( ((F-Fbar)**2).mean())
            X = this_sim.slope_rights
            Fbar.shape=Fbar.size
            ax.errorbar( X, Fbar, yerr=std)


    outname='plots_to_sort/spectra_%s_%s.png'%(plotname,field)
    fig.savefig(outname)
    plt.close(fig)
    print(outname)






