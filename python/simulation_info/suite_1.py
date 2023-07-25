from GL import *
import simulation
reload(simulation)
import sim_colors

if 1:
    sim_ms = nar(['half','1','2','3','4','5','6'])
    sim_ms_f = nar([0.5,1,2,3,4,5,6])
    sim_ma = nar(['half','1','2'])
    sim_ma_f = nar([0.5,1,2])
if 0:
    #for kludging
    sim_ms = nar(['half','1','2','3'])
    sim_ms_f = nar([0.5,1,2,3])
    sim_ma = nar(['half','1','2'])
    sim_ma_f = nar([0.5,1,2])

#color_by_mach = {'half':'c','1':'m','2':'b','3':'g','5':'r'}
color_by_mach = {'half':'red','1':'orange','2':'g','3':'b','4':'violet','5':'brown','6':'black'}
line_by_alf_mach  = {'half':':','1':'--','2':'-'}
marker_by_alf_mach = {'half':'.','1':'^','2':'s'}

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
        color[sim]=color_by_mach[ms]
        linestyle[sim]=line_by_alf_mach[ma]
        marker[sim] = marker_by_alf_mach[ma]
        #glyph = color[sim]+linestyle[sim]
        tdyn[sim] = 0.5/sim_ms_f[nmach]

#auto gen, don't touch
simlist = nar([ '%s_%s'%(ms,ma) for ms in sim_ms for ma in sim_ma])

markerlist = nar([ marker['%s_%s'%(ms,ma)] for ms in sim_ms for ma in sim_ma])
colorlist  = nar([ color['%s_%s'%(ms,ma)] for ms in sim_ms for ma in sim_ma])
linelist  = nar([ linestyle['%s_%s'%(ms,ma)] for ms in sim_ms for ma in sim_ma])

def vals_from_sim(sim):
    ms,ma = sim.split("_")
    if ms == 'half':
        ms = 0.5
    if ma == 'half':
        ma = 0.5
    ms=float(ms)
    ma=float(ma)
    return ms,ma

ms_list=[]
ma_list=[]
for sim in simlist:
    ms,ma = vals_from_sim(sim)
    ms_list.append( ms)
    ma_list.append(ma)
ms_list=nar(ms_list)
ma_list=nar(ma_list)
Ms = dict(zip(simlist,ms_list))
Ma = dict(zip(simlist,ma_list))

def lrange(*args):
    return list(range(*args))
three_half_range = list( range(70,85))+list(range(86,93))
analysis_frames={
    "half_half":lrange(11,30)+lrange(31,32),"half_1":lrange(11,30)+lrange(31,32),"half_2":lrange(11,30)+lrange(31,32),
    "1_half":lrange(11,30)+lrange(31,32),"1_1":lrange(11,30)+lrange(31,32),"1_2":lrange(11,30)+lrange(31,32),
    "2_half":lrange(65,84)+lrange(85,86),"2_1":lrange(11,30)+lrange(31,32),"2_2":lrange(11,26)+lrange(27,32),
    #"3_half":lrange(72,93),"3_1":lrange(56,75),"3_2":lrange(20,40),
    "3_half":lrange(70,85)+lrange(90,30)+lrange(91,93),"3_1":lrange(53,74)+lrange(75,77),"3_2":lrange(9,39)+lrange(40,41),
    #'4_half':lrange(12,19), '4_1':lrange(12,22),'4_2':lrange(12,25),
    #'4_half':lrange(15,45), '4_1':lrange(15,52),'4_2':lrange(15,52),
    '4_half':lrange(14,31)+lrange(32,45), '4_1':lrange(16,23)+lrange(24,52),'4_2':lrange(15,37)+lrange(38,52),
    "5_half":lrange(11,40)+lrange(41,60),"5_1":lrange(4,37)+lrange(38,49),"5_2":lrange(4,46)+lrange(48,60),"5_3":lrange(5,59),
    '6_half':lrange(17,39)+lrange(40,46), '6_1':lrange(15,30)+lrange(31,52),'6_2':lrange(16,30)+lrange(31,52)}

for sim in simlist:
    simulation.sim(sim, data_location=dl.sim_dir_base+sim, product_location=dl.product_dir_base+sim, ms=Ms[sim], ma=Ma[sim],
                   color=color[sim],linestyle=linestyle[sim],marker=marker[sim],
                   framelist=analysis_frames[sim])
