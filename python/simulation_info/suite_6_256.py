from GL import *
import simulation
reload(simulation)
import sim_colors

if 1:
    sim_ms = nar(['1','2','3','4','5','6','7','8'])
    sim_ms_f = sim_ms.astype('float')
    sim_ma = nar(['0'])
    sim_ma_f = sim_ma.astype('float')



#auto gen, don't touch
simlist = nar([ '%s_%s'%(ms,ma) for ms in sim_ms for ma in sim_ma])

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
    longname[sim] = 'a%02d_Ms%0.1f_Ma%0.1f_256'%(counter, ms, ma)
    longnamelist.append(longname[sim])
    longname_from_key['a%02d'%counter]=longname[sim]
ms_list=nar(ms_list)
ma_list=nar(ma_list)
Ms = dict(zip(simlist,ms_list))
Ma = dict(zip(simlist,ma_list))


def launch_script():
    for sim in simlist:
        print('turb_maker.py -n %s -s %0.1f -a %0.1f -d 256'%(longname[sim][:3], Ms[sim], Ma[sim]))

analysis_frames={}
for sim in simlist:
    analysis_frames[sim] = range(1,11)

data_base = "/anvil/scratch/x-ux454321/Paper83/suite_6_256_athena/"
product_base = "/anvil/scratch/x-ux454321/Paper83/suite_6_256_athena/products"

for sim in simlist:
    simulation.sim(longname[sim], data_location="%s/%s"%(data_base,longname[sim]), 
                   product_location="%s/%s"%(product_base,longname[sim]), ms=Ms[sim], ma=Ma[sim],
                   color='k',linestyle='-',marker='*',
                   all_frames=analysis_frames[sim],framelist=analysis_frames[sim], code='Athena')
