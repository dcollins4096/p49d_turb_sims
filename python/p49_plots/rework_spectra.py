
from GL import *


from GL import *

import simulation


def rework_spectra(simlist):


    Nzones, kspace = dt.dpy('nzones.h5',['Nzones','kspace'])
    for ns,sim in enumerate(simlist):
        this_sim = simulation.corral[sim]
        for frame in this_sim.all_frames:
            k3d, density = dt.dpy('%s/DD%04d.products/power_density.h5'%(this_sim.product_location,frame), ['k','power'])
            k3d, magnetic  = dt.dpy('%s/DD%04d.products/power_magnetic.h5'%(this_sim.product_location,frame), ['k','power'])
            k3d, velocity = dt.dpy('%s/DD%04d.products/power_velocity.h5'%(this_sim.product_location,frame), ['k','power'])
            avg_dens_name='%s/DD%04d.products/avg_power_density.h5'%(this_sim.product_location,frame)
            avg_magn_name='%s/DD%04d.products/avg_power_magnetic.h5'%(this_sim.product_location,frame)
            avg_velo_name='%s/DD%04d.products/avg_power_velocity.h5'%(this_sim.product_location,frame)
            names=[avg_dens_name,avg_velo_name,avg_magn_name]
            specs=[density,velocity,magnetic]
            for n in [0,1,2]:
                print(names[n])
                fptr=h5py.File(names[n],'w')
                fptr['k']=kspace
                fptr['power']=specs[n]/Nzones
                fptr.close()
