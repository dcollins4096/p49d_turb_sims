from starter1 import *
import yt
import downsample.volavg as volavg
import filament.hessian as hessian
import tools.pcolormesh_helper as pch
import tools.davetools as dt
import filament.color_games as cg
reload(cg)

def asinh(arr,sigma):
    const=0.25
    c0=const*np.arcsinh(np.abs(arr)/(const*sigma))
    c0/=c0.max()
    return c0
def pproj(estuff,outname):
    fig,axes=plt.subplots(4,3,figsize=(12,12))
    scale_proj=estuff.proj
    scale_proj-=scale_proj.min()
    scale_proj/=scale_proj.max()
    print(scale_proj.min(),scale_proj.max())
    axes[0][0].imshow(cg.color_packer(scale_proj, scale_proj,scale_proj))
    bins=np.linspace(-100,100,1000)
    kwargs={'histtype':'step','bins':bins}
    axes[0][1].hist( estuff.eP0f,label='eP0',color='m',**kwargs)
    axes[0][1].hist( estuff.eP1f,label='eP1',color='c',**kwargs)
    pch.simple_phase( estuff.eP0f, estuff.eP1f, ax=axes[0][2])


    sigma = np.std(estuff.eP0)
    const=0.25
    c0=asinh(estuff.eP0,sigma)
    c1=asinh(estuff.eP1,sigma)
    zeros=np.zeros_like(c1)
    axes[1][0].imshow(cg.color_packer(c0, zeros, zeros))
    axes[1][1].imshow(cg.color_packer(zeros, c1, zeros))
    axes[1][2].imshow(cg.color_packer(c0, c1, zeros))

    paxis=0
    print(estuff.e0.shape)
    pe0=estuff.e0.sum(axis=paxis)
    pe1=estuff.e1.sum(axis=paxis)
    pe2=estuff.e2.sum(axis=paxis)
    sigma = np.std(estuff.e1)
    const=0.25
    ee0=asinh(pe0,sigma)
    ee1=asinh(pe1,sigma)
    ee2=asinh(pe2,sigma)
    axes[2][0].imshow(cg.color_packer(ee0, zeros, zeros))
    axes[2][1].imshow(cg.color_packer(zeros, ee2, zeros))
    axes[2][2].imshow(cg.color_packer(ee0, ee2, zeros))

    axes[0][1].hist( pe0.flatten(),label='e0',color='r',**kwargs)
    axes[0][1].hist( pe1.flatten(),label='e1',color='g',**kwargs)
    axes[0][1].hist( pe2.flatten(),label='e2',color='b',**kwargs)

    xx1,yy1=estuff.eP0f.flatten(), pe0.flatten()
    xx2,yy2=estuff.eP1f.flatten(), pe2.flatten()
    pfit1=np.polyfit(xx1,yy1,1)
    pfit2=np.polyfit(xx2,yy2,1)
    print("fit1",pfit1)
    print("fit2",pfit2)
    pch.simple_phase( xx1,yy1, ax=axes[3][0])
    pch.simple_phase( xx2,yy2, ax=axes[3][1])
    xi = nar([xx1.min(),xx1.max()])
    xi2 = nar([xx2.min(),xx2.max()])
    axes[3][0].plot( xi,pfit1[0]*xi+pfit1[1])
    axes[3][1].plot( xi2,pfit2[0]*xi2+pfit2[1])
    #axes[3][0].plot( xi, xi-33)
    #axes[3][0].plot([-100,20],[-100,20])
    #axes[3][1].plot([-25,100],[-25,100])
    #axes[2][2].contour(Q,levels=[15])

    if 1:
        x1 = cg.color_vec(c0, [1.0,0.5,0.0])
        x2 = cg.color_vec(ee0, [0.0,0.5,1.0])
        axes[3][2].imshow( x1+x2)
    if 0:
        Q = np.abs(estuff.eP0)+np.abs(estuff.eP1)
        Qs=asinh(Q,sigma)
        axes[3][2].imshow( cg.color_packer(Qs,Qs,Qs))
        axes[0][1].hist(Q.flatten(), label='Q',color='k',**kwargs)

    axes[0][1].legend(loc=0)
    fig.savefig(outname)

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

def get_cubes(sim,frame):
    ds = yt.load("/data/cb1/Projects/P49_EE_BB/%s/DD%04d/data%04d"%(sim, frame, frame))
    print('get cg')
    cg = ds.covering_grid(0, [0.0]*3, [512]*3)
    dds = cg.dds
    rho_full = cg["density"].v
    rho = volavg.volavg(rho_full,rank=3,refine_by=2)
    return rho_full, rho

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
    def doproj(self,axis):
        self.rho_periodic=hessian.extend_periodic(self.arr,1)
        self.proj = self.rho_periodic.sum(axis=axis)
        self.hess_proj = hessian.hessian(self.proj,self.dds)
        self.evals_proj,self.evecs_proj=hessian.eigensystem(self.hess_proj)
        hessian.sort_hessian(self.evals_proj,self.evecs_proj)

        self.eP0= self.evals_proj[...,0]
        self.eP1= self.evals_proj[...,1]
        #self.eP2= self.evals_proj[...,2]
        self.eP0f=self.evals_proj[...,0].flatten()
        self.eP1f=self.evals_proj[...,1].flatten()
        #self.eP2f=self.evals_proj[...,2].flatten()
        print('hoot')

        self.done=True

