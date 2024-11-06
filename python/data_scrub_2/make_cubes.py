
from GL import *
import yt
import sim_colors
import downsample.downsampler as DOON
import simulation
import data_locations as dl
reload(DOON)

import simulation_info.all_sims as all_sims

def make_cubes(sim_list):
    for sim in sim_list:
        this_sim=simulation.corral[sim]
        #this_sim.load()
        if this_sim.B_nom > 0:
            do_magnetic=True
        else:
            do_magnetic=False
        for frame in this_sim.ann_frames:
            print('Load',sim,frame)

            #print(source_fname)
            already_got_this_dir=False
            dir_128='/anvil/scratch/x-ux454321/p83_turbulence/Athena/Cubes/128/%s'%sim
            if not os.path.exists( dir_128 ):
                os.makedirs( dir_128 )
            frame_dir = "%s/DD%04d"%(dir_128,frame)
            if not os.path.exists( frame_dir ):
                os.makedirs( frame_dir )


            #dest_fname = "%s/%s/DD%04d/data%04d.cube.h5"%(sim_colors.cloudbreak_128, sim, frame, frame)
            refine_by = 4
            dest_fname = "%s/cube%04d"%(frame_dir,frame)
            files = glob.glob("%s*"%dest_fname)
            if len(files) > 0:
                print("Got some files", files)
                continue
            ds = this_sim.load_ds(frame)
            DOON.downsample_and_write(ds,dest_fname, write_hdf5=False,write_fits=True, refine_by=refine_by, do_magnetic = do_magnetic)


