
# Started and abandoned.  Showed that the vorticity is zero and the driving is as expected.
#got as far as showing its zero for one sim.
#Pretty sure we're solid.


from  GL import *
import spectra_tools as st

#import simulation_info.suite_liltest
#reload(simulation_info.suite_liltest)
#import simulation_info.suite_1
#reload(simulation_info.suite_1)
import simulation
import simulation_info.all_sims as all_sims

import make_vorticity_spectra as mvs
reload(mvs)
this_list = all_sims.lists['suite1']
this_list = ['6_2']


if 0:
    mvs.make_spec(this_list,frames=-1)

for sim in this_list:
    this_sim=simulation.corral[sim]
    Nzones, kspace = dt.dpy('nzones.h5',['Nzones','kspace'])
    frames = this_sim.all_frames[-1:]
    for frame in frames:
        oober = st.short_oober(this_sim.data_location, frame=frame, product_directory=this_sim.product_location, 
                               simname=sim)

        if 'vx' not in dir():
            print('read and avg')
            ax = oober.fft(frame,'x-acceleration')
            print('read and avg')
            ay = oober.fft(frame,'y-acceleration')
            print('read and avg')
            az = oober.fft(frame,'z-acceleration')
            print('read and avg')
            powerx,nzones,ff = st.shell_average_only(ax*np.conj(ax))
            print('read and avg')
            powery,nzones,ff = st.shell_average_only(ax*np.conj(ay))
            print('read and avg')
            powerz,nzones,ff = st.shell_average_only(ax*np.conj(az))

    kx,ky,kz = ff._kk
    if 'Omega' not in dir():
        Omega = ax*(ky*az-kz*ay)+ ay*(kz*ax-kx*az)+ az*(kx*ay-ky*ax)
        print('average')
        power_omega,nzones_o,ff_o = st.shell_average_only(Omega*np.conj(Omega))
    fig,axes=plt.subplots(1,1)
    ax0=axes

    ax0.plot(kspace, powerx)
    ax0.plot(kspace, power_omega)
    ax0.set(xscale='log',yscale='log')
    fig.savefig('%s/vorticity_test_%s'%(plotdir,oober.simname))
