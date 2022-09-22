from GL import *
import h5py
import get_all_quantities as gaq
from cycler import cycler
from queb3 import powerlaw_fit as plfit
reload(gaq)
#Makes quan plots for each frame and each sim
#quan=gaq.all_quan_from_outputlog('/scratch/00369/tg456484/Paper49d_moresims/zd01_M10_MA1_512_quan')
#'/scratch/00369/tg456484/Paper49d_moresims/zd01_M10_MA1_512_quan'
# out_prefix = 'zd01'
generaldir='/data/cb1/Projects/P49_EE_BB/plots/general'
simdic=["half_half","half_1","half_2","1_half","1_1","1_2","2_half","2_1","2_2","3_half","3_1","3_2"]
frames=[range(12,31),range(12,31),range(12,31),range(12,31),range(12,31),range(12,31),range(66,85),range(12,31),range(12,31),range(72,91),range(57,75),range(21,40)]
alfmachavg=np.zeros(12)
sonicmachavg=np.zeros(12)
brmsavgsim=np.zeros(12)
v2bavgsim=np.zeros(12)
masought=[0.5,1.0,2.0,0.5,1.0,2.0,0.5,1.0,2.0,0.5,1.0,2.0]
mssought=[0.5,0.5,0.5,1.0,1.0,1.0,2.0,2.0,2.0,3.0,3.0,3.0]
for i in range(12):
  direc='/data/cb1/Projects/P49_EE_BB/frbs/%s'%simdic[i]
  plotdir='/home/dccollins/PigPen'
#    all_files =gaq.all_quan_from_files(direc,frames[i])
#    files_to_use = all_files #all_files[slice(None,None,30)]+all_files[-1:]
#    quan = gaq.return_average_quantities(files_to_use)
  v2avg=[]
  b2avg=[]
  vaavg=[]
  maavg=[]
  timeavg=[]
  brmsavg=[]
  v2bavg=[]
  

#  afm=0
#  sma=0


  for frame in frames[i]:
    quan=h5py.File('%s/data%04d.AverageQuantities.h5'%(direc,frame))

    #These are code quantities, so you do not want a sqrt(4pi) if comparing to other code values
    v2 = np.sqrt(quan['vx_std'][:]**2+quan['vy_std'][:]**2+quan['vz_std'][:]**2)
    b2 = np.sqrt(quan['bx_std'][:]**2+quan['by_std'][:]**2+quan['bz_std'][:]**2+
                 quan['bx_avg'][:]**2+quan['by_avg'][:]**2+quan['bz_avg'][:]**2)
    va = np.sqrt(quan['alf_x_std'][:]**2+quan['alf_y_std'][:]**2+quan['alf_z_std'][:]**2)
    bstd = np.sqrt(quan['bx_std'][:]**2+quan['by_std'][:]**2+quan['bz_std'][:]**2)
    bt = np.sqrt(quan['bx_avg'][:]**2+quan['by_avg'][:]**2+quan['bz_avg'][:]**2)
    #Ma should be v2/btot (v/B is how we defined ma in code initializer) (btot b2 without the bstd)
    print(v2)
    ma=v2/bt
    vrms2=v2**2
    brms=bstd
    v2b=vrms2/bt
    
    v2avg=np.append(v2avg,v2)
    b2avg=np.append(b2avg,b2)
    vaavg=np.append(vaavg,va)
    maavg=np.append(maavg,ma)
    v2bavg=np.append(v2bavg,v2b)
    brmsavg=np.append(brmsavg,brms)
    timeavg=np.append(timeavg,quan['time'][:])

    v2btime=np.mean(v2bavg)
    v2bavgsim[i]=v2btime
    brmstime=np.mean(brmsavg)
    brmsavgsim[i]=brmstime
#    matime=np.mean(ma)
#    mstime=np.mean(v2)
#    afm+=matime
#    print('afm=')
#    print(afm)
#    sma+=mstime
#    afm2+=matime2
#   plt.clf()
#   plt.plot(quan['time'][:],quan['vx_avg'][:],marker="*",label='vx_avg')
#   plt.plot(quan['time'][:],quan['vy_avg'][:],marker="*",label='vy_avg')
#   plt.plot(quan['time'][:],quan['vz_avg'][:],marker="*",label='vz_avg')
#   plt.xlabel('t');plt.ylabel(r'$\langle v_i \rangle$')
#   plt.legend(loc=0)
#   plt.savefig('%s/%s_time_vi_avg_%04d.png'%(plotdir,simdic[i],frame))

#   plt.clf()
#   plt.plot(quan['time'][:],quan['vx_std'][:],marker="*",label='vx_rms')
#   plt.plot(quan['time'][:],quan['vy_std'][:],marker="*",label='vy_rms')
#   plt.plot(quan['time'][:],quan['vz_std'][:],marker="*",label='vz_rms')
#   plt.xlabel('t');plt.ylabel(r'$\sqrt{ \langle v_i^2\rangle - \langle v_i \rangle^2}$')
#   plt.legend(loc=0)
#   plt.savefig('%s/%s_time_vi_rms_%04d.png'%(plotdir,simdic[i],frame))

#   plt.clf()
#   plt.plot(quan['time'][:],quan['bx_avg'][:],marker="*",label='bx_avg')
#   plt.plot(quan['time'][:],quan['by_avg'][:],marker="*",label='by_avg')
#   plt.plot(quan['time'][:],quan['bz_avg'][:],marker="*",label='bz_avg')
#   plt.xlabel('t');plt.ylabel(r'$\langle B_i \rangle$')
#   plt.legend(loc=0)
#   plt.savefig('%s/%s_time_bi_avg_%04d.png'%(plotdir,simdic[i],frame))

#   plt.clf()
#   plt.plot(quan['time'][:],quan['bx_std'][:],marker="*",label='bx_rms')
#   plt.plot(quan['time'][:],quan['by_std'][:],marker="*",label='by_rms')
#   plt.plot(quan['time'][:],quan['bz_std'][:],marker="*",label='bz_rms')
#   plt.xlabel('t');plt.ylabel(r'$\sqrt{ \langle B_i^2\rangle - \langle B_i \rangle^2}$')
#   plt.legend(loc=0)
#   plt.savefig('%s/%s_time_bi_rms_%04d.png'%(plotdir,simdic[i],frame))
#   plt.clf()
#  afm/=len(frames[i])
#  sma/=len(frames[i])
#  afm2/=len(frames[i])
#  alfmachavg[i]=afm
#  sonicmachavg[i]=sma
  alfmachavg[i]=np.mean(maavg)
  sonicmachavg[i]=np.mean(v2avg)
  print('alfmachavg =')
  print(alfmachavg)
  print('sonicmachavg =')
  print(sonicmachavg)
  plt.plot(timeavg,vaavg,marker="*")
  plt.xlabel('t');plt.ylabel(r'$Va_{\rm{rms}}$')
  plt.savefig('%s/%s_time_va.png'%(plotdir,simdic[i]))
  plt.clf()
  plt.plot(timeavg,v2avg,marker="*")
  plt.xlabel('t');plt.ylabel(r'$V_{\rm{rms}}$')
  plt.savefig('%s/%s_time_mach.png'%(plotdir,simdic[i]))
  plt.clf()
  plt.plot(timeavg,b2avg,marker="*")
  plt.xlabel('t');plt.ylabel(r'$\sqrt{\langle B_i B_i\rangle}$')
  plt.savefig('%s/%s_time_b.png'%(plotdir,simdic[i]))
  plt.clf()
  plt.plot(timeavg,maavg,marker="*")
  plt.xlabel('t');plt.ylabel(r'$Ma_{\rm{rms}}$')
  plt.savefig('%s/%s_time_ma.png'%(plotdir,simdic[i]))
  plt.clf()
  plt.plot(brmstime,v2btime,marker="*")
  plt.xlabel(r'$b_{\rm{rms}}$');plt.ylabel(r'${v_{\rm{rms}}**2}/{b_{\rm{tot}}}$')
  plt.savefig('%s/%s_brms_v2b.png'%(plotdir,simdic[i]))
  plt.clf()
  plt.plot(timeavg,vaavg,label=r'$V_a$')
  plt.plot(timeavg,v2avg,label=r'$V_{rms}$')
  plt.plot(timeavg,b2avg,label=r'$B_{rms}$')
  plt.plot(timeavg,maavg,label=r'$M_a$')
  plt.xlabel('t')
  plt.legend(loc=0)
  plt.savefig('%s/%s_all_quan.png'%(plotdir,simdic[i]))
  plt.clf()
  
msmax=np.max(sonicmachavg)+1.0
mamax=np.max(alfmachavg)+1.0
xms=np.linspace(0,msmax,20)
xma=np.linspace(0,mamax,20)
cm=plt.get_cmap('gist_rainbow')
fig = plt.figure()
axis=fig.add_subplot(111)
#axis.set_prop_cycle(color=[cm(1.*i/len(labels)) for i in range(len(labels))]) #set colormap to cycle colors i
axis.set_prop_cycle(cycler(color=['r','g','b','r','g','b','r','g','b','r','g','b'])+cycler(marker=['.','.','.','x','x','x','*','*','*','1','1','1']))
for i in range(len(simdic)):
    axis.plot(mssought[i], sonicmachavg[i],ms=5,label='%s'%simdic[i])
plt.xlabel(r'$M_{s}$ Sought')
plt.ylabel(r'$M_{s}$ Measured')
axis.plot(xms, xms,label='Ideal 1 to 1', c='k',marker='')
plt.legend(loc=0)
plt.savefig('%s/Ms_paramspace.png'%(plotdir))
plt.clf()

fig2, ax2 = plt.subplots(1,1)
cm=plt.get_cmap('gist_rainbow')
fig = plt.figure()
axis=fig.add_subplot(111)
#axis.set_prop_cycle(color=[cm(1.*i/len(labels)) for i in range(len(labels))]) #set colormap to cycle colors i
axis.set_prop_cycle(cycler(color=['r','g','b','r','g','b','r','g','b','r','g','b'])+cycler(marker=['.','.','.','x','x','x','*','*','*','1','1','1']))
for i in range(len(simdic)):
    axis.plot(masought[i], alfmachavg[i],ms=5,label='%s'%simdic[i])
    ax2.plot(sonicmachavg[i], alfmachavg[i])
fig2.savefig('%s/msmameas.png'%plotdir)
plt.xlabel(r'$M_{a}$ Sought')
plt.ylabel(r'$M_{a}$ Measured')
axis.plot(xma, xma,label='Ideal 1 to 1', c='k',marker='')
plt.legend(loc=0)
plt.savefig('%s/Ma_v_div_va_paramspace.png'%(plotdir))
plt.clf()

#cm=plt.get_cmap('gist_rainbow')
#fig = plt.figure()
#axis=fig.add_subplot(111)
#axis.set_prop_cycle(cycler(color=['r','g','b','r','g','b','r','g','b','r','g','b'])+cycler(marker=['.','.','.','x','x','x','*','*','*','1','1','1']))
#for i in range(len(simdic)):
#    axis.plot(masought[i], alfmachavg2[i],ms=5,label='%s'%simdic[i])
#plt.xlabel(r'$M_{a}$ Sought')
#plt.ylabel(r'$M_{a}$ Measured')
#axis.plot(xma, xma,label='Ideal 1 to 1', c='k',marker='')
#plt.legend(loc=0)
#plt.savefig('%s/Ma_v_div_B_paramspace.png'%(generaldir))
#plt.clf()

cm=plt.get_cmap('gist_rainbow')
fig = plt.figure()
axis=fig.add_subplot(111)
#fitrange=np.zeros(2)
#fitrange[0]=brmsavgsim[0]
#fitrange[1]=brmsavgsim[-1]
#slopes=list(plfit(brmsavgsim,v2bavgsim,fitrange))  #returns [slope,amp,res]
bmax=np.max(brmsavgsim)+.25
xslope=np.linspace(0,bmax,10)
m,b=np.polyfit(brmsavgsim,v2bavgsim,1)
yslope=(m*xslope)+b
print('slope is')
print(m)
axis.plot(xslope,yslope,label=r'$slope$ %0.2f'%m,c='k',marker='')
axis.set_prop_cycle(cycler(color=['r','g','b','r','g','b','r','g','b','r','g','b'])+cycler(marker=['.','.','.','x','x','x','*','*','*','1','1','1']))
for i in range(len(simdic)):
    axis.plot(brmsavgsim[i], v2bavgsim[i],ms=5,label='%s'%simdic[i])
plt.xlabel(r'$B_{rms}$')
plt.ylabel(r'$V_{a}^{2}/B_{tot}$')
plt.legend(loc=0)
plt.savefig('%s/brms_v2b.png'%(plotdir))
plt.clf()

Fptr = h5py.File("ma_ms_dict.h5",'w')
Fptr.create_dataset("Ma", data=alfmachavg)
Fptr.create_dataset("Ms", data=sonicmachavg)
Fptr.close()
