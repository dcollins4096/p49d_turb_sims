from starter1 import *

plt.close('all')


import filament.fake_filament as fil
import filament.hessian as hessian
import filament.tools as htools
reload(fil)
reload(htools)

print('good morning')
plotdir = "%s/PigPen"%os.environ['HOME']

rho = fil.make_random_1(3, 32)

if 'esystem' not in dir() or True:
    esystem=htools.eigen_stuff(rho, cut_periodic=True)

if not esystem.done:
    esystem.do()

if 0:
    esystem.hist_evals("%s/ev2"%plotdir)
if 1:
    esystem.proj("%s/proj"%plotdir)


if 0:
    fig,axes=plt.subplots(1,3, figsize=(8,4))
    axes[0].imshow(rho.sum(axis=0),origin='lower',interpolation='nearest')
    axes[1].imshow(rho.sum(axis=1),origin='lower',interpolation='nearest')
    axes[2].imshow(rho.sum(axis=2),origin='lower',interpolation='nearest')
    fig.tight_layout()
    fig.savefig('%s/fake1_filaments'%plotdir)

if 0:
    if 'logrho_periodic' not in dir():
        print('extend')
        #logrho_periodic = hessian.extend_periodic(np.log(rho),1)
        logrho_periodic = hessian.extend_periodic(np.log(rho),1)
    if 'hess' not in dir():
        print('hessian')
        dds=[1,1,1]
        hess = hessian.hessian(logrho_periodic,dds)
    if 'evals' not in dir():
        print('eigen')
        evals,evecs=hessian.eigensystem(hess)
        hessian.sort_hessian(evals,evecs)

    if 'e0' not in dir():
        e0=evals[:,:,:,0]
        e1=evals[:,:,:,1]
        e2=evals[:,:,:,2]

    if 1:
        fig,axes=plt.subplots(1,3)
        print('histograms')
        #bins=np.linspace(-1e5,1e5,128)
        bins=None
        axes[0].hist(e0.flatten(),bins=bins, label='L0', histtype='step')
        axes[0].hist(e1.flatten(),bins=bins, label='L1', histtype='step')
        axes[0].hist(e2.flatten(),bins=bins, label='L2', histtype='step')
        axes[0].legend(loc=0)
        axes[0].set(yscale='log')
        fig.savefig('%s/derp'%plotdir)
