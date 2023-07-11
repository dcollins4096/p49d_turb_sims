
from GL import *
import queb3
import get_all_quantities as gaq
#import sim_colors
#
# Read restricted set.
#
import simulation as sim
from collections import defaultdict
if 'quan3' not in dir():
    quan3={}
    quan_time={}
def read(sim_name, clobber=False):
    if sim_name in quan3 and not clobber:
        print("Not reading twice", sim_name)
        return

    this_sim = sim.corral[sim_name]
    quan3[sim_name]={}
    quan_time[sim_name]={}
    v2avg=[]
    msavg=[]
    maavg=[]
    for frame in this_sim.framelist:
        fname = '%s/DD%04d.products/data%04d.AverageQuantities.h5'%(this_sim.product_location,frame,frame)
        if not os.path.exists(fname):
            print("missing",fname)
            continue
        quan4=h5py.File(fname,'r')
        try:
            for field in quan4:
                if field in quan_time[sim_name]:
                    quan_time[sim_name][field] = np.concatenate([quan_time[sim_name][field],quan4[field][()]])
                else:
                    quan_time[sim_name][field] = quan4[field][()]
            v2 = np.sqrt(quan4['vx_std'][:]**2+quan4['vy_std'][:]**2+quan4['vz_std'][:]**2)
            bt = np.sqrt(quan4['bx_avg'][:]**2+quan4['by_avg'][:]**2+quan4['bz_avg'][:]**2)
            ma=v2/bt
            v2avg=np.append(v2avg,v2)
            maavg=np.append(maavg,ma)
        except:
            raise
        finally:
            quan4.close()
    quan3[sim_name]['maavg'] =np.mean(maavg)
    quan3[sim_name]['msavg'] =np.mean(v2avg)

