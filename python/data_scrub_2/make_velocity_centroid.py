
from GL import *
import queb3
reload(queb3)
import simulation

def make_all(simlist):
    print('simlist',simlist)
    for sim in simlist:
        print(sim)
        this_sim = simulation.corral[sim]

        prefix = this_sim.name
        pack = queb3.simulation_package( directory=this_sim.data_location,frames=this_sim.all_frames,prefix=prefix, 
                                        product_directory=this_sim.product_location, simname=sim, code=this_sim.code)
        for frame in pack.frames:
            pack.make_frbs(frame, do_magnetic=False,do_velocity=True)

def image_centroid(simlist):
    for sim in simlist:
        this_sim = simulation.corral[sim]
        for frame in this_sim.ann_frames:
            fig,ax=plt.subplots(2,3)
            for nax,axis in enumerate('xyz'):
                cen = "%s/DD%04d.products/DD%04d_velocity_centroid_%s.fits"%(this_sim.product_location, frame, frame, axis)
                var = "%s/DD%04d.products/DD%04d_velocity_variance_%s.fits"%(this_sim.product_location, frame, frame, axis)
                if os.path.exists(cen):
                    d1 = pyfits.open(cen)[0].data
                    ax[0][nax].imshow(d1)
                    ax[1][nax].hist(d1.flatten())
            outname = '%s/plots/%s_centroid_%04d.png'%(os.environ['HOME'],sim,frame)
            fig.savefig(outname)
            print(outname)
            plt.close(fig)


