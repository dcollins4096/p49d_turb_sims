
from GL import *
import downsample.volavg as volavg
import filament.hessian as hessian

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

if 1:
    fig,axes=plt.subplots(2,2)
    print('histograms')
    bins=np.linspace(-1e5,1e5,128)
    e0=evals[:,:,:,0].flatten()
    e1=evals[:,:,:,1].flatten()
    e2=evals[:,:,:,2].flatten()
    axes[0][0].hist(e0 ,bins=bins, label='L0', histtype='step')
    axes[0][0].hist(e1 ,bins=bins, label='L1', histtype='step')
    axes[0][0].hist(e2 ,bins=bins, label='L2', histtype='step')
    axes[0][0].legend(loc=0)
    axes[0][0].set(yscale='log')
    fig.savefig( "%s/evals"%plotdir)

if 0:
    fig,axes=plt.subplots(2,2)
    print('der')
    bins=np.linspace(-1e3,1e3,128)
    for ne,eee in enumerate([e0,e1,e2]):
        es = np.abs(eee)
        print('sort')
        es.sort()
        es=es[::10]
        print( "%0.2e"%es[y<0.75].max())
        y = np.arange(es.size)/es.size
        axes[0][0].plot(es,y,label='L%d'%ne)
    axes[0][0].set(xscale='log',yscale='linear')
    axes[0][0].axhline(0.05)
    axes[0][0].axhline(0.95)
    print('save')
    fig.savefig('%s/cuml'%plotdir)
