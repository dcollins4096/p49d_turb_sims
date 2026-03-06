from GL import *
import simulation
reload(simulation)
import sim_colors

if 1:
    sim_ms = nar(['1','1.5','2','2.5','3','3.5','4','4.5','5','5.5','6','6.5','7','7.5','8','8.5'])
    sim_ms_f = sim_ms.astype('float')
    sim_ma = nar(['0'])
    sim_ma_f = sim_ma.astype('float')


#color_by_mach = {'half':'c','1':'m','2':'b','3':'g','5':'r'}
color_by_mach = {'1':'cyan','2':'red','3':'orange','4':'g','5':'b','6':'violet','7':'brown','8':'black'}
line_by_alf_mach  = {'0.5':':','1.5':'--','3':'-', '0':':-'}
marker_by_alf_mach = {'0.5':'.','1.5':'^','3':'s', '0':'*'}

plot_order=[]
color={}
linestyle={}
marker={}
glyph={}
tdyn={}
for nma,ma in enumerate(sim_ma):
    for nmach,ms in enumerate(sim_ms):
        sim="%s_%s"%(ms,ma)
        plot_order.append(sim)
        color[sim]=color_by_mach.get(ms, 'k')
        linestyle[sim]=line_by_alf_mach.get(ma,'-')
        marker[sim] = marker_by_alf_mach.get(ma,'*')
        #glyph = color[sim]+linestyle[sim]
        tdyn[sim] = 0.5/sim_ms_f[nmach]

#auto gen, don't touch
simlist = nar([ '%s_%s'%(ms,ma) for ms in sim_ms for ma in sim_ma])

markerlist = nar([ marker['%s_%s'%(ms,ma)] for ms in sim_ms for ma in sim_ma])
colorlist  = nar([ color['%s_%s'%(ms,ma)] for ms in sim_ms for ma in sim_ma])
linelist  = nar([ linestyle['%s_%s'%(ms,ma)] for ms in sim_ms for ma in sim_ma])

def vals_from_sim(sim):
    ms,ma = sim.split("_")
    ms=float(ms)
    ma=float(ma)
    return ms,ma

ms_list=[]
ma_list=[]
longname = {}
long_from_key={}
list_from_key={}
long_simlist = []
sim_from_key={}
for counter,sim in enumerate(simlist):
    ms,ma = vals_from_sim(sim)
    key = 'd%02d'%counter
    ms_list.append( ms)
    ma_list.append(ma)
    longname[sim] = '%s_Ms%0.1f_Ma%0.1f_256'%(key, ms, ma)
    long_simlist.append(longname[sim])
    list_from_key[key]=[longname[sim]]
    long_from_key[key] = longname[sim]
    sim_from_key[key]=sim
ms_list=nar(ms_list)
ma_list=nar(ma_list)
Ms = dict(zip(simlist,ms_list))
Ma = dict(zip(simlist,ma_list))

def launch_script():
    for sim in simlist:
        print('turb_maker.py -n %s -s %0.1f -a %0.1f -d 256'%(longname[sim][:3], Ms[sim], Ma[sim]))

analysis_frames={}
for sim in simlist:
    analysis_frames[sim] = list(range(20,1002))
    #analysis_frames[sim] = [1,30]
#analysis_frames[sim_from_key['c28']] = range(1,19)
data_location = "/anvil/scratch/x-ux454321/Paper83/suite_5_256_long"
product_location = "/anvil/scratch/x-ux454321/Paper83/suite_5_products"
for sim in simlist:
    simulation.sim(longname[sim], 
                   data_location="%s/%s"%(data_location,longname[sim]), 
                   product_location="%s/%s"%(product_location,longname[sim]), ms=Ms[sim], ma=Ma[sim],
                   color=color[sim],linestyle=linestyle[sim],marker=marker[sim],
                   all_frames=analysis_frames[sim],framelist=analysis_frames[sim], code='Enzo')
