from GL import *

import simulation
import time as tt
#https://github.com/idies/pyJHTDB/blob/master/examples/isotropic_spectra_3D.ipynb
def spectra_test(sim,frame):

    print('read')
    this_sim=simulation.corral[sim]
    vxhat=dt.dpy("%s/DD%04d.products/fft_velocity_x.float32"%(this_sim.product_location,frame),['velocity_x'])[0]
    vyhat=dt.dpy("%s/DD%04d.products/fft_velocity_y.float32"%(this_sim.product_location,frame),['velocity_y'])[0]
    vzhat=dt.dpy("%s/DD%04d.products/fft_velocity_z.float32"%(this_sim.product_location,frame),['velocity_z'])[0]
    dim=vxhat.shape[0]//2

    uu_fft=(np.abs(vxhat[::2,::2,::2])/dim**3)**2
    vv_fft=(np.abs(vyhat[::2,::2,::2])/dim**3)**2
    ww_fft=(np.abs(vzhat[::2,::2,::2])/dim**3)**2
    L=2*np.pi


    k_end=int(dim/2)
    rx=np.array(range(dim))-dim/2+1
    rx=np.roll(rx,int(dim/2)+1)

    start = tt.time()
    r=np.zeros((rx.shape[0],rx.shape[0],rx.shape[0]))
    print('horse around')
    for i in range(rx.shape[0]):
        for j in range(rx.shape[0]):
            r[i,j,:]=rx[i]**2+rx[j]**2+rx[:]**2
    r=np.sqrt(r)
    print(tt.time() - start)

    dx=2*np.pi/L
    k=(np.array(range(k_end))+1)*dx

    start = tt.time()
    bins=np.zeros((k.shape[0]+1))
    for N in range(k_end):
        if N==0:
            bins[N]=0
        else:
            bins[N]=(k[N]+k[N-1])/2    
    bins[-1]=k[-1]

    inds = np.digitize(r*dx, bins, right=True)
    spectrum=np.zeros((k.shape[0]))
    bin_counter=np.zeros((k.shape[0]))

    for N in range(k_end):
        spectrum[N]=np.sum(uu_fft[inds==N+1])+np.sum(vv_fft[inds==N+1])+np.sum(ww_fft[inds==N+1])
        bin_counter[N]=np.count_nonzero(inds==N+1)

    spectrum=spectrum*2*np.pi*(k**2)/(bin_counter*dx**3)
    print(tt.time() - start)
    fig,axes= plt.subplots(1,2)
    ax=axes[0]
    ax1=axes[1]
    ax.plot(k,spectrum)
    s2=3e-1*k**(-5/3)
    s2*=spectrum[1]/s2[1]
    ax.plot(k,s2)
    ax.set_xscale('log')
    ax.set_yscale('log')


    ax1.plot( k, 2*np.pi*(k**2)/(bin_counter*dx**3))
    fig.savefig('plots_to_sort/test')

spectra_test("4_2",40)

