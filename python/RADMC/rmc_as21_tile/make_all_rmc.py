"""
Make many queb data products.

*  queb3.simulation_package points to a simulation.  
   BoxSize is in units of 64 zones.
*  
      
"""
from GL import *
from davetools import *
import queb3
reload(queb3)
frames=[1]
sim_dir = "rmc_as21_tile"
plot_dir= os.environ['HOME']+'/PigPen/'

pack = queb3.simulation_package( directory=sim_dir,frames=frames,prefix="as21_window",
                                dataset_name='as21_tilewindow_DD',frbname="./",plotdir=plot_dir)
pack_tile = queb3.simulation_package( directory=sim_dir,frames=frames,prefix="as21_window_tile",
                                dataset_name='as21_tilewindow_DD',frbname="./",plotdir=plot_dir)
for frame in frames:
    theta_phi =  [[0,0,'zish']]
    for theta,phi,lab in theta_phi:
        print(theta,phi,lab)
        this_proj=pack.read_queb(frame=frame,theta=theta,phi=phi)
        this_proj_tile=pack_tile.read_queb(frame=frame,theta=theta,phi=phi)

        this_proj_tile.tile(2)
        this_proj_tile.compute_harmonic_products()
        this_proj_tile.untile(2)

        this_proj.compute_harmonic_products()

        pack.image_fields_6way(frame,axis="0",ts=this_proj,
                               field_list=['E','B', 'Q','U','T'],theta=theta,phi=phi)
        pack_tile.image_fields_6way(frame,axis="0",ts=this_proj_tile,
                               field_list=['E','B', 'Q','U','T'],theta=theta,phi=phi)

