from GL import *
import davetools as dt

sim_ms = nar(['half','1','2','3','4','5','6'])
sim_ms_f = nar([0.5,1,2,3,4,5,6])
#sim_ms = nar(['half','1','2','3']); print('kludge; no sim 5')
sim_ma = nar(['half','1','2'])
sim_ma_f = nar([0.5,1,2])

#auto gen, don't touch
simlist = nar([ '%s_%s'%(ms,ma) for ms in sim_ms for ma in sim_ma])


rm = dt.rainbow_map(4)
#color_by_mach = {'half':'c','1':'m','2':'b','3':'g','5':'r'}
color_by_mach = {'half':'red','1':'orange','2':'g','3':'b','4':'violet','5':'brown','6':'black'}
line_by_alf_mach  = {'half':':','1':'--','2':'-'}
marker_by_alf_mach = {'half':'.','1':'^','2':'s'}


#color=dict(zip(simlist,color_list))
#linestyle = dict(zip(simlist,line_list))
#marker = dict(zip(simlist,marker_list))
#frames = dict(zip(simlist,framelist))
#glyph = dict(zip(simlist,glyph_list))

#color_list=nar(['r','g','b','r','g','b','r','g','b','r','g','b','r','g','b'])
#line_list=nar(['-','-','-','-.','-.','-.','--','--','--',':',':',':', '-','-','-'])
#glyph_list = nar([a+b for a,b in zip(color_list,line_list)])
#marker_list = nar(['.','.','.','s','s','s','^','^','^','*','*','*','o','o','o'])

three_half_range = list( range(70,85))+list(range(86,93))
def lrange(*args):
    return list(range(*args))
framedict={
    "half_half":lrange(11,32),"half_1":lrange(11,32),"half_2":lrange(11,32),
    "1_half":lrange(11,32),"1_1":lrange(11,32),"1_2":lrange(11,32),
    "2_half":lrange(65,86),"2_1":lrange(11,32),"2_2":lrange(11,32),
    #"3_half":lrange(72,93),"3_1":lrange(56,75),"3_2":lrange(20,40),
    "3_half":three_half_range,"3_1":lrange(53,77),"3_2":lrange(9,41),
    "5_half":lrange(3,37),"5_1":lrange(4,28),"5_2":lrange(4,46),"5_3":lrange(5,59),
    #'4_half':lrange(12,19), '4_1':lrange(12,22),'4_2':lrange(12,25),
    '4_half':lrange(15,45), '4_1':lrange(15,52),'4_2':lrange(15,52),
    '6_half':lrange(12,19), '6_1':lrange(12,22),'6_2':lrange(12,25)}
frames=framedict
#framedict['5_half']=range(11,13); print('kludge in framedict')



#
# General setup.
#
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

markerlist = nar([ marker['%s_%s'%(ms,ma)] for ms in sim_ms for ma in sim_ma])
colorlist  = nar([ color['%s_%s'%(ms,ma)] for ms in sim_ms for ma in sim_ma])
linelist  = nar([ linestyle['%s_%s'%(ms,ma)] for ms in sim_ms for ma in sim_ma])



framelist=[framedict[sim] for sim in simlist]

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


