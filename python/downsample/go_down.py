from GL import *
import yt
import sim_colors
import downsampler as DOON
reload(DOON)

counter=0
for sim in sim_colors.simlist:
    for frame in sim_colors.frames[sim]:
        counter += 1
        source_fname = "%s/%s/DD%04d/data%04d"%(sim_colors.cloudbreak_base, sim, frame, frame)
        if not os.path.exists(source_fname):
            continue
        dir_128 =  "%s/%s/DD%04d/"%(sim_colors.cloudbreak_128, sim, frame)
        print(source_fname)
        already_got_this_dir=False
        if not os.path.exists( dir_128 ):
            os.mkdir( dir_128 )
        else:
            already_got_this_dir=True

        #Maybe we want do not repeat ourselves.
        if already_got_this_dir:
            continue

        dest_fname = "%s/%s/DD%04d/data%04d.cube.h5"%(sim_colors.cloudbreak_128, sim, frame, frame)
        ds = yt.load(source_fname)
        DOON.downsample_and_write(ds,dest_fname)
        print(ds)


