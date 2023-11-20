from starter1 import *
import yt
from dtools import volavg
import fourier_tools_py3.fourier_filter as Filter
import brunt_tools as bt
reload(bt)
plt.close('all')
plotdir='%s/PigPen/'%(os.environ['HOME'])




#sim = '4_1'
#frame = 35
#sim='1_1'
#frame=31
sim='6_2'
frame=51
#sim = ms_ma ms = 123456 ma = half,1,2

#get cubes; rho_full is straight off disk.  rho is downsampled by 2
if 'rho' not in dir():
    print('load and cg')
    rho_full, rho = bt.get_cubes(sim,frame)
    print(sim,frame)

if 'ftool' not in dir():
    ftool = bt.fft_tool(rho)
    ftool.do3()
    ftool.do2(projax=0)

if 1:
    bt.plot_fft(ftool,"%s/ffts_%s_%04d"%(plotdir,sim,frame))

if 1:
    outname="%s/fft_full_%s_n%04d"%(plotdir,sim,frame)
    bt.plot_brunt(ftool, outname=outname,method='full')

    outname="%s/fft_norm_%s_n%04d"%(plotdir,sim,frame)
    bt.plot_brunt(ftool, outname=outname,method='norm')

    outname="%s/fft_range_%s_n%04d"%(plotdir,sim,frame)
    bt.plot_brunt(ftool, outname=outname,method='range', fitrange=[4,25])

if 0:
    #try brunt.  Works to a factor of 2, not bad.
    sigma_rho = (ftool.rho**2).sum().real
    sigma_fft = (ftool.power_1d3).sum().real
    #sigma_fft = (ftool.power).sum().real/ftool.power.size**2
    print("rho %0.2e fft %0.2e ratio %0.2e"%(sigma_rho,sigma_fft,(sigma_fft)/sigma_rho))

    sigma_col = ((ftool.rho2)**2).sum().real
    sigma_2d_fft = ftool.power_1d2.sum().real
    print("col %0.2e fft %0.2e ratio %0.2e"%(sigma_col, sigma_2d_fft, sigma_col/sigma_2d_fft))
    Rinv = ( ftool.k2d*ftool.power_1d2)[1:].sum()/(ftool.power_1d2)[1:].sum()
    Rinv_actual = ftool.power_1d3.sum()/ftool.power_1d2.sum()
    sigma_Brunt = sigma_col*Rinv
    sigma_Brunt_actual = sigma_col*Rinv_actual
    print( sigma_Brunt/sigma_rho)
    print( sigma_Brunt_actual/sigma_rho)
    fig,ax=plt.subplots(1,1)
    Q=ftool.k2d*ftool.power_1d2
    P=ftool.power_1d3
    ax.plot( ftool.k2d, Q)
    ax.plot( ftool.k3d, P)
    ax.set(xscale='log',yscale='log')
    fig.savefig("%s/test"%plotdir)
