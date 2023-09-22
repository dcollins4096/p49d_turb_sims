from GL import *

import simulation

def plot(simlist, Nzones2d=1, Nzones3d=1):

    fig,ax0=plt.subplots(1,1)
    #fig,axes=plt.subplots(2,2)
    #ax0=axes[0][0];ax1=axes[0][1]
    #ax2=axes[1][0];ax3=axes[1][1]
    all_gamma=[]
    
    for ns,sim in enumerate(simlist):
        this_sim = simulation.corral[sim]
        this_sim.load()
        package=this_sim.return_queb3_package()
        this_proj=package.read_queb(30,'y')
        print(this_proj.T)

        xvals2 = this_sim.avg_spectra['k2d']/(2*np.pi)
        xvals3 = this_sim.avg_spectra['k3d']
        TT= this_sim.avg_spectra['ClTTy'] /Nzones2d[1:]
        rho= this_sim.avg_spectra['density']/Nzones3d
        v = TT[20]/rho[20]
        rho=rho*v
        #pdb.set_trace()
        print(v/np.pi)
        #fit_range =this_sim.get_fitrange(xvals)
        #mask = (xvals > fit_range[0])*(xvals < fit_range[1])

        kwargs={'color':this_sim.color,'linestyle':this_sim.linestyle}
        ax0.plot( xvals2, TT, c='r')
        ax0.plot( xvals3, rho, c='b')



    ax0.set(xscale='log',yscale='log')


    fig.savefig('plots_to_sort/spectra_test')
