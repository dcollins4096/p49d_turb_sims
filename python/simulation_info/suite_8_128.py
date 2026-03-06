from GL import *
import simulation
reload(simulation)
import sim_colors

if 1:
    sim_ms_f = np.arange(0.5,5.25,0.25)
    sim_ma_f = nar([0])



#auto gen, don't touch
simlist = nar([ '%0.2f_%0.2f'%(ms,ma) for ms in sim_ms_f for ma in sim_ma_f])

def vals_from_sim(sim):
    ms,ma = sim.split("_")
    ms=float(ms)
    ma=float(ma)
    return ms,ma

ms_list=[]
ma_list=[]
longname = {}
longnamelist = []
longname_from_key={}
for counter,sim in enumerate(simlist):
    ms,ma = vals_from_sim(sim)
    ms_list.append( ms)
    ma_list.append(ma)
    longname[sim] = 'b%02d_Ms%0.1f_Ma%0.1f_128'%(counter, ms, ma)
    longnamelist.append(longname[sim])
    longname_from_key['b%02d'%counter]=longname[sim]
ms_list=nar(ms_list)
ma_list=nar(ma_list)
Ms = dict(zip(simlist,ms_list))
Ma = dict(zip(simlist,ma_list))


def launch_script():
    for sim in simlist:
        print('turb_maker.py -n %s -s %0.1f -a %0.1f -d 256'%(longname[sim][:3], Ms[sim], Ma[sim]))

analysis_frames={}
for sim in simlist:
    analysis_frames[sim] = range(1000)

data_base = "/anvil/scratch/x-ux454321/Paper83/suite_8_128_athena/"
product_base = "/anvil/scratch/x-ux454321/Paper83/suite_8_128_athena/products"

for sim in simlist:
    simulation.sim(longname[sim], data_location="%s/%s"%(data_base,longname[sim]), 
                   product_location="%s/%s"%(product_base,longname[sim]), ms=Ms[sim], ma=Ma[sim],
                   color='k',linestyle='-',marker='*',
                   all_frames=analysis_frames[sim],framelist=analysis_frames[sim], code='Athena')
