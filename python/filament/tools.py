from starter1 import *

import filament.hessian as hessian
import tools.pcolormesh_helper as pch
import tools.davetools as dt

def hist_evals(estuff, outname,bins=None):
    fig,axes=plt.subplots(1,1)
    print('histograms')
    if bins==True:
        bins=np.linspace(-1e5,1e5,128)
    axes.hist(estuff.e0f ,bins=bins, label='L0', histtype='step')
    axes.hist(estuff.e1f ,bins=bins, label='L1', histtype='step')
    axes.hist(estuff.e2f ,bins=bins, label='L2', histtype='step')
    axes.legend(loc=0)
    axes.set(yscale='log')
    fig.savefig(outname)
    plt.close(fig)
    print("Made",outname)

def proj(estuff,outname=None):
    fig,axes=plt.subplots(4,3,figsize=(12,12))
    for axis in [0,1,2]:
        xyz='xyz'[axis]
        p=axes[0][axis].imshow( estuff.arr.sum(axis=axis))
        axes[0][axis].set(title='rho '+xyz)
        fig.colorbar(p,ax=axes[0][axis])
        p=axes[1][axis].imshow( estuff.e0.sum(axis=axis))
        axes[1][axis].set(title='e0 '+xyz)
        fig.colorbar(p,ax=axes[1][axis])
        p=axes[2][axis].imshow( estuff.e1.sum(axis=axis))
        axes[2][axis].set(title='e1 '+xyz)
        fig.colorbar(p,ax=axes[2][axis])
        p=axes[3][axis].imshow( estuff.e2.sum(axis=axis))
        fig.colorbar(p,ax=axes[3][axis])
        axes[3][axis].set(title='e2 '+xyz)

    fig.tight_layout()
    fig.savefig(outname)

def contour(estuff,outname=None):
    def doer(arr,axis):
        sl = [slice(None),slice(None),slice(None)]
        sl[axis]=12
        return arr[tuple(sl)]
    fig,ax=plt.subplots(1,1,figsize=(8,8))
    kwargs={'origin':'lower','interpolation':'nearest'}
    axis=0
    xyz='xyz'[axis]
    p=ax.imshow( doer(estuff.arr,axis), **kwargs)
    ax.set(title='rho '+xyz)
    fig.colorbar(p,ax=ax)
    ax.contour( doer(estuff.e0,axis),levels=[-0.25])
    fig.savefig(outname)

def color_game(estuff,outname=None):
    def doer(arr,axis):
        sl = [slice(None),slice(None),slice(None)]
        sl[axis]=12
        b= arr[tuple(sl)]
        c = np.log(np.abs(b))
        c /= c.max()-c.min()
        c -= c.min()
        return c

    fig,ax=plt.subplots(2,3,figsize=(8,8))
    ext=dt.extents()
    ext(estuff.e0)
    ext(estuff.e1)
    ext(estuff.e2)
    a = np.array(ext.minmax)
    edge = np.max(np.abs(a))
    bins = np.linspace(-edge,edge,128)

    def packer(a,b,c):
        return np.stack([a,b,c]).transpose()
    rho= doer(estuff.arr,0)
    a1 = doer(estuff.e0,0)
    a2 = doer(estuff.e1,0)
    a3 = doer(estuff.e2,0)
    zero=np.zeros_like(a1)
    #ax.imshow(y)
    ax[0][0].imshow( packer(a1,zero,zero))
    ax[0][1].imshow( packer(zero,a2,zero))
    ax[0][2].imshow( packer(zero,zero,a3))
    ax[1][0].imshow( packer(a1,a2,a3))
    ax[1][2].imshow( packer(a1,zero,a3))
    #pch.simple_phase(estuff.e0.flatten(),estuff.e1.flatten(),bins=[bins,bins], ax=axes[0])
    #pch.simple_phase(estuff.e1.flatten(),estuff.e2.flatten(),bins=[bins,bins], ax=axes[1])
    #pch.simple_phase(estuff.e0.flatten(),estuff.e2.flatten(),bins=[bins,bins], ax=axes[2])
    #axes[0].set(title='E1 vs E0')
    #axes[1].set(title='E2 vs E1')
    #axes[2].set(title='E2 vs E0')
    fig.tight_layout()
    fig.savefig(outname)
    plt.close(fig)
    fig,ax=plt.subplots(1,1,figsize=(8,8))
    ax.imshow(packer(a1,zero,zero))
    fig.savefig('%s_e0'%outname)
    ax.clear()
    ax.imshow(packer(a2,zero,zero))
    fig.savefig('%s_e1'%outname)
    ax.clear()
    ax.imshow(packer(a3,zero,zero))
    fig.savefig('%s_e2'%outname)
    plt.close(fig)
    def packer2(a,vec):
        return np.stack([a*vec[0],a*vec[1],a*vec[2]]).transpose()
    fig,ax=plt.subplots(1,3,figsize=(8,4))
    v1=packer2(a1,[1,0.5,0])
    ax[0].imshow(v1)
    v2=packer2(a3,[0,0.5,1])
    ax[1].imshow(v2)
    ax[2].imshow(v1+v2)
    #ax[1][0].imshow(packer(rho,rho,rho))
    #ax[1][0].imshow(packer(rho,rho,rho))
    #ax[1][1].imshow(packer(a1,a3,rho))
    fig.savefig('%s_game2'%outname)
    plt.close(fig)
    fig,ax=plt.subplots(2,2,figsize=(8,8))
    ax[0][0].imshow(packer(a1,zero,zero))
    ax[0][1].imshow(packer(zero,a3,zero))
    ax[1][0].imshow(packer(zero,zero,rho))
    ax[1][1].imshow(packer(a1,a3,rho))
    fig.savefig('%s_game3'%outname)
    plt.close(fig)





def correlations(estuff,outname=None):

    fig,axes=plt.subplots(1,3,figsize=(8,4))
    ext=dt.extents()
    ext(estuff.e0)
    ext(estuff.e1)
    ext(estuff.e2)
    a = np.array(ext.minmax)
    edge = np.max(np.abs(a))
    bins = np.linspace(-edge,edge,128)

    pch.simple_phase(estuff.e0.flatten(),estuff.e1.flatten(),bins=[bins,bins], ax=axes[0])
    pch.simple_phase(estuff.e1.flatten(),estuff.e2.flatten(),bins=[bins,bins], ax=axes[1])
    pch.simple_phase(estuff.e0.flatten(),estuff.e2.flatten(),bins=[bins,bins], ax=axes[2])
    axes[0].set(title='E1 vs E0')
    axes[1].set(title='E2 vs E1')
    axes[2].set(title='E2 vs E0')
    fig.tight_layout()
    fig.savefig(outname)



def doslice(estuff,outname=None):
    fig,axes=plt.subplots(4,3,figsize=(12,12))
    def doer(arr,axis):
        return arr.sum(arr,axis=axis)
    def doer(arr,axis):
        sl = [slice(None),slice(None),slice(None)]
        sl[axis]=12
        return arr[tuple(sl)]
    for axis in [0,1,2]:
        xyz='xyz'[axis]
        kwargs={'origin':'lower','interpolation':'nearest'}
        p=axes[0][axis].imshow( doer(estuff.arr,axis), **kwargs)
        axes[0][axis].set(title='rho '+xyz)
        fig.colorbar(p,ax=axes[0][axis])
        p=axes[1][axis].imshow( doer(estuff.e0,axis),**kwargs)
        axes[1][axis].set(title='e0 '+xyz)
        fig.colorbar(p,ax=axes[1][axis])
        p=axes[2][axis].imshow( doer(estuff.e1,axis),**kwargs)
        axes[2][axis].set(title='e1 '+xyz)
        fig.colorbar(p,ax=axes[2][axis])
        p=axes[3][axis].imshow( doer(estuff.e2,axis),**kwargs)
        fig.colorbar(p,ax=axes[3][axis])
        axes[3][axis].set(title='e2 '+xyz)

    fig.tight_layout()
    fig.savefig(outname)
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

