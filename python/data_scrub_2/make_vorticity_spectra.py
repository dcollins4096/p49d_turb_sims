#
# Started and abandoned.  Showed that the vorticity is zero and the driving is as expected.
#
from GL import *
from queb3 import powerlaw_fit as plfit
import spectra_tools as st

import simulation

import fourier_tools_py3.fourier_filter as Filter

def vorticity(oober,frame):

    print(oober.simname,frame)
    vx = oober.fft(frame,'velocity_x')
    vy = oober.fft(frame,'velocity_y')
    vz = oober.fft(frame,'velocity_z')
    powerx,nzones,ff = st.shell_average_only(vx*np.conj(vx))
    powery,nzones,ff = st.shell_average_only(vx*np.conj(vy))
    powerz,nzones,ff = st.shell_average_only(vx*np.conj(vz))
    Nzones, kspace = dt.dpy('nzones.h5',['Nzones','kspace'])

    fig,ax=plt.subplots(1,1)
    ax0=ax

    ax0.plot(kspace, powerx)
    ax0.set(xscale='log',yscale='log')
    fig.savefig('%s/vorticity_test_%s'%(plotdir,oober.simname))


def make_spec(simlist, frames=None):
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
    for nsim,sim in enumerate(simlist):
        this_sim=simulation.corral[sim]


        if frames is None:
            frames = this_sim.all_frames
        elif frames == -1:
            frames = this_sim.all_frames[-1:]
        Nzones, kspace = dt.dpy('nzones.h5',['Nzones','kspace'])
        for frame in frames:
            oober = st.short_oober(this_sim.data_location, frame=frame, product_directory=this_sim.product_location, 
                                   simname=sim)
            if 1:
                st.MakeAccelSpectra(oober,frame)
                k3d, accel = dt.dpy('%s/DD%04d.products/power_acceleration.h5'%(this_sim.product_location,frame), ['k','power'])
                avg_accel_name='%s/DD%04d.products/avg_power_acceleration.h5'%(this_sim.product_location,frame)
                fptr=h5py.File(avg_accel_name,'w')
                fptr['k']=kspace
                fptr['power']=accel/Nzones
                fptr.close()
            if 0:
                vorticity(oober,frame)
