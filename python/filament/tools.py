from starter1 import *

import filament.hessian as hessian
class eigen_stuff():
    def __init__(self,array,dds=None, cut_periodic=False):
        self.dds=dds
        if dds==None:
            self.dds = np.ones(len(array.shape))

        self.arr=array
        self.done=False
        self.cut_periodic=cut_periodic
    def do(self):
        self.rho_periodic=hessian.extend_periodic(self.arr,1)
        self.hess = hessian.hessian(self.rho_periodic,self.dds)
        self.evals,self.evecs=hessian.eigensystem(self.hess)
        hessian.sort_hessian(self.evals,self.evecs)
        if self.cut_periodic:
            self.evals = self.evals[1:-1,1:-1,1:-1,:]

        self.e0=self.evals[:,:,:,0]
        self.e1=self.evals[:,:,:,1]
        self.e2=self.evals[:,:,:,2]
        self.e0f=self.evals[:,:,:,0].flatten()
        self.e1f=self.evals[:,:,:,1].flatten()
        self.e2f=self.evals[:,:,:,2].flatten()
        self.done=True
    def hist_evals(self, outname,bins=None):
        fig,axes=plt.subplots(1,1)
        print('histograms')
        if bins==True:
            bins=np.linspace(-1e5,1e5,128)
        axes.hist(self.e0f ,bins=bins, label='L0', histtype='step')
        axes.hist(self.e1f ,bins=bins, label='L1', histtype='step')
        axes.hist(self.e2f ,bins=bins, label='L2', histtype='step')
        axes.legend(loc=0)
        axes.set(yscale='log')
        fig.savefig(outname)
        plt.close(fig)
        print("Made",outname)
    def proj(self,outname=None):
        fig,axes=plt.subplots(4,3,figsize=(12,12))
        for axis in [0,1,2]:
            xyz='xyz'[axis]
            p=axes[0][axis].imshow( self.arr.sum(axis=axis))
            axes[0][axis].set(title='rho '+xyz)
            fig.colorbar(p,ax=axes[0][axis])
            p=axes[1][axis].imshow( self.e0.sum(axis=axis))
            axes[1][axis].set(title='e0 '+xyz)
            fig.colorbar(p,ax=axes[1][axis])
            p=axes[2][axis].imshow( self.e1.sum(axis=axis))
            axes[2][axis].set(title='e1 '+xyz)
            fig.colorbar(p,ax=axes[2][axis])
            p=axes[3][axis].imshow( self.e2.sum(axis=axis))
            fig.colorbar(p,ax=axes[3][axis])
            axes[3][axis].set(title='e2 '+xyz)

        fig.tight_layout()
        fig.savefig(outname)

