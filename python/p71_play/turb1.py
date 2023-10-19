

#
# Make hessian on downsampled cube.
#


from starter1 import *
import downsample.volavg as volavg
import filament.hessian as hessian
import tools.pcolormesh_helper as pch
import yt

plotdir = "%s/PigPen"%(os.environ['HOME'])


sim = '4_1'
frame = 35



if 'ds' not in dir():
    print('load and cg')
    ds = yt.load("/data/cb1/Projects/P49_EE_BB/%s/DD%04d/data%04d"%(sim, frame, frame))
    cg = ds.covering_grid(0, [0.0]*3, [512]*3)
    dds = cg.dds
    rho_full = cg["density"].v
    rho = volavg.volavg(rho_full,rank=3,refine_by=2)
    rho256=rho
if 'logrho_periodic' not in dir():
    print('extend')
    logrho_periodic = hessian.extend_periodic(np.log(rho),1)
if 'hess' not in dir():
    print('hessian')
    hess = hessian.hessian(logrho_periodic,dds)
if 'evals' not in dir():
    print('eigen')
    evals,evecs=hessian.eigensystem(hess)
    hessian.sort_hessian(evals,evecs)

if 'e0' not in dir():
    e0=evals[:,:,:,0].flatten()
    e1=evals[:,:,:,1].flatten()
    e2=evals[:,:,:,2].flatten()

if 0:
    fig,axes=plt.subplots(1,3)
    print('histograms')
    bins=np.linspace(-1e5,1e5,128)
    axes[0].hist(e0 ,bins=bins, label='L0', histtype='step')
    axes[0].hist(e1 ,bins=bins, label='L1', histtype='step')
    axes[0].hist(e2 ,bins=bins, label='L2', histtype='step')
    axes[0].legend(loc=0)
    axes[0].set(yscale='log')

    axes[1].hist(e0[e0>0] ,bins=bins, label='L0', histtype='step')
    axes[1].hist(e1[e1>0] ,bins=bins, label='L1', histtype='step')
    axes[1].hist(e2[e2>0] ,bins=bins, label='L2', histtype='step')
    axes[1].legend(loc=0)
    axes[1].set(yscale='log', xscale='log', title='e>0')
    axes[2].hist(-e0[e0<0] ,bins=bins, label='L0', histtype='step')
    axes[2].hist(-e1[e1<0] ,bins=bins, label='L1', histtype='step')
    axes[2].hist(-e2[e2<0] ,bins=bins, label='L2', histtype='step')
    axes[2].legend(loc=0)
    axes[2].set(yscale='log', xscale='log', title='e<0')

    fig.tight_layout()
    fig.savefig( "%s/evals"%plotdir)

if 0:
    fig,axes=plt.subplots(1,3)
    print('phase')
    bins=np.linspace(-1e5,1e5,128)
    pch.simple_phase(e0,e1,bins=[bins,bins], ax=axes[0])
    pch.simple_phase(e1,e2,bins=[bins,bins], ax=axes[1])
    pch.simple_phase(e0,e2,bins=[bins,bins], ax=axes[2])
    axes[0].set(title='E1 vs E0')
    axes[1].set(title='E2 vs E1')
    axes[2].set(title='E2 vs E0')
    fig.tight_layout()
    fig.savefig("%s/phase"%plotdir)

if 1:
    fig,axes=plt.subplots(1,3, figsize=(8,4))
    print('play')
    bins=np.linspace(-1e5,1e5,128)
    pch.simple_phase(e0,e1-e0,bins=[bins,bins], ax=axes[0])
    pch.simple_phase(e2,e2-e0,bins=[bins,bins], ax=axes[1])
    pch.simple_phase(e1,e0+e1+e2,bins=[bins,bins], ax=axes[2])
    axes[0].set(title='E1 vs E0')
    axes[1].set(title='E2 vs E1')
    axes[2].set(title='E2 vs E0')
    fig.tight_layout()
    fig.savefig("%s/play"%plotdir)







