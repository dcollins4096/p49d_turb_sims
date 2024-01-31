from GL import *
import yt
import downsampler as DOON
reload(DOON)
import simulation
import simulation_info.all_sims as all_sims

sim_list=all_sims.lists['suite1']
for sim in sim_list:
    this_sim = simulation.corral[sim]
    for frame in this_sim.ann_frames[-1:]:
        base_dir="/data/cb1/Share/MHDTurbulence"
        dir_128 =  "%s/%s/"%(base_dir, sim)
        already_got_this_dir=False
        if not os.path.exists( dir_128 ):
            os.mkdir( dir_128 )
        else:
            already_got_this_dir=True

        #Maybe we want do not repeat ourselves.
        #if already_got_this_dir:
        #    continue

        dest_fname = "%s/data%04d.cube.h5"%(dir_128,frame)
        if os.path.exists(dest_fname):
            continue

        ds = this_sim.load_ds(frame)
        #DOON.downsample_and_write(ds,dest_fname, refine_by=1, write_hdf5=True,write_fits=False)
        cg = ds.covering_grid(0,[0.0]*3,ds.domain_dimensions)
        fptr = h5py.File(dest_fname,"w")
        for in_field_name in ['Density',
                              'x-velocity','y-velocity','z-velocity',
                              'Bx','By','Bz']: #['Density','x-velocity','y-velocity','z-velocity','BxF','ByF','BzF']:
            print("Read %s"%in_field_name)
            field_to_store = cg[in_field_name]
            print("Write %s"%in_field_name)
            fptr.create_dataset(in_field_name, data=field_to_store)
        fptr.close()


