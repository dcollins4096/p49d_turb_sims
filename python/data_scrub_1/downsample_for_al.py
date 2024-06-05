
from GL import *
import yt
import sim_colors
import downsample.downsampler as DOON
import simulation
import data_locations as dl
reload(DOON)

import simulation_info.all_sims as all_sims
sim_list = all_sims.lists['suite1']
for sim in sim_list:
    this_sim=simulation.corral[sim]
    this_sim.load()
    for frame in this_sim.ann_frames[-1:]:
        print('Load',sim,frame)
        source_fname = "%s/DD%04d/data%04d"%(this_sim.data_location, frame, frame)
        if not os.path.exists(source_fname):
            print(source_fname)
            print('missing')
            continue
        if 0:
            dir_128 =  "%s/%s/DD%04d/"%(dl.p58_dir, sim, frame)
            refine_by = 2
        else:
            dir_128 =  "%s/512/%s/DD%04d/"%(dl.p58_dir, sim, frame)
            refine_by = 1

        #print(source_fname)
        already_got_this_dir=False
        if not os.path.exists( dir_128 ):
            os.mkdir( dir_128 )
        else:
            already_got_this_dir=True

        #Maybe we want do not repeat ourselves.
        if already_got_this_dir:
            continue

        #dest_fname = "%s/%s/DD%04d/data%04d.cube.h5"%(sim_colors.cloudbreak_128, sim, frame, frame)
        dest_fname = "%s/cube_%04d"%(dir_128,frame)
        ds = yt.load(source_fname)
        DOON.downsample_and_write(ds,dest_fname, write_hdf5=False,write_fits=True, refine_by=refine_by)


