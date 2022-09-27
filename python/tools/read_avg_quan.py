
from GL import *
import queb3
import get_all_quantities as gaq
import sim_colors
if 0:
    if 'quand' not in dir():
        quand = {}
    for i, sim in enumerate(sim_colors.simlist):
        if sim in quand:
            continue
        print("=== %s ==="%sim)
        sim_dir = "/data/cb1/Projects/P49_EE_BB/%s"%sim
        ol = "%s/OutputLog"%(sim_dir)
        quand[sim] = gaq.get_quantities_and_rms(ol)
#
# Read restricted set.
#

quan3={}
for ns, sim in enumerate(sim_colors.simlist):
    quan3[sim]={}
    v2avg=[]
    msavg=[]
    maavg=[]
    for frame in sim_colors.framelist[ns]:
        fname = '%s/Products/%s/DD%04d.products/data%04d.AverageQuantities.h5'%(sim_colors.cloudbreak_base,sim,frame,frame)
        if not os.path.exists(fname):
            print("missing",fname)
            continue
        quan4=h5py.File(fname,'r')
        v2 = np.sqrt(quan4['vx_std'][:]**2+quan4['vy_std'][:]**2+quan4['vz_std'][:]**2)
        bt = np.sqrt(quan4['bx_avg'][:]**2+quan4['by_avg'][:]**2+quan4['bz_avg'][:]**2)
        ma=v2/bt
        v2avg=np.append(v2avg,v2)
        maavg=np.append(maavg,ma)
    quan3[sim]['maavg'] =np.mean(maavg)
    quan3[sim]['msavg'] =np.mean(v2avg)
