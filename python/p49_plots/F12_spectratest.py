from GL import *

import simulation

if 'stuff' not in dir():
    stuff={}
import fourier_tools_py3.fourier_filter as Filter
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
        frame=30
        this_proj=package.read_queb(frame,'y')
        this_proj.compute_harmonic_products()
        k_k=this_sim.avg_spectra['k2d']
        rho2=this_proj.T
        fft2 = np.fft.fftn( rho2 )
        power2=fft2*np.conjugate(fft2)
        ff2 = Filter.FourierFilter(power2)
        power_1d2 = np.array([power2[ff2.get_shell(bin)].sum() for bin in range(ff2.nx)])
        power_1d2 /= rho2.size
        Nzones2 = np.array([ff2.get_shell(bin).sum() for bin in range(ff2.nx)])
        kspace2a=ff2.get_shell_k()
        dk=kspace2a[1]
        kspace2=2*np.pi*(kspace2a+0.5*dk)

        #print(k_k[:3])
        #print(2*np.pi*(kspace2[:3]+0.5/512))
        #print(2*np.pi/k_k[:3])
        k3d = 2*np.pi*(this_sim.avg_spectra['k3d']+0.5*dk)
        spec_3d = this_sim.all_spectra[frame]['density']

        if 1:
            if os.path.exists('nzones.h5') :
                Nzones, kspace = dt.dpy('nzones.h5',['Nzones','kspace'])
            else:
                power=np.zeros([512,512,512])
                ff = Filter.FourierFilter(power)
                print('count')
                Nzones = np.array([ff.get_shell(bin).sum() for bin in range(ff.nx)])
                kspace=ff.get_shell_k()
                kspace=2*np.pi*(kspace+0.5/512)
                fptr=h5py.File('nzones.h5','w')
                fptr['Nzones']=Nzones
                fptr['kspace']=kspace
                fptr.close()

                print('done')
        if 0:
            ds = this_sim.load_ds(frame)
            Nz = 512
            cg=ds.covering_grid(0,[0.0]*3,[Nz]*3)
            rho=cg['density'].v
            fft1 = np.fft.fftn( rho )
            power=fft1*np.conjugate(fft1)
            ff = Filter.FourierFilter(power)
            power_1d = np.array([power[ff.get_shell(bin)].sum() for bin in range(ff.nx)])
            Nzones1 = np.array([ff.get_shell(bin).sum() for bin in range(ff.nx)])
            kspace=ff.get_shell_k()
            kspace=2*np.pi*(kspace+0.5/512)


        fig,ax=plt.subplots(1,1)
        ax.plot( k_k, this_proj.ClTT,c='r')
        ax.plot( kspace2, power_1d2/Nzones2,c='g')
        ax.plot( k3d, 512**2*spec_3d/Nzones, c='b')
        #ax.plot( kspace, power_1d/Nzones1, c='k')
        ax.set(xscale='log',yscale='log')
        fig.savefig('plots_to_sort/ft_test')
        print('wut')

