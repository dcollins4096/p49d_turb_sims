from dtools.starter1 import *
import yt
from downsample import volavg
import fourier_tools_py3.fourier_filter as Filter
import brunt_tools as bt
reload(bt)
plt.close('all')

#DEFINE PLOTDIR
#plotdir='%s/PigPen/'%(os.environ['HOME'])

sim='6_2'
frame=51

#get cubes; rho_full is straight off disk.  rho is downsampled by 2
if 'rho' not in dir():
    print('load and cg')
    rho_full, rho = bt.get_cubes(sim,frame)
    print(sim,frame)

if 'ftool' not in dir() or True:
    #X1,Y1,Z1=10,20,30
    #rho = rho_full[X1:X1+256,Y1:Y1+256,Z1:Z1+256]
    ftool = bt.fft_tool(rho)
    ftool.do3()
    ftool.do2(projax=0)

if 1:
    outname="%s/fft_full_%s_n%04d"%(plot_dir,sim,frame)
    bt.plot_brunt(ftool, outname=outname,method='full')


if 0:
    plotdir="./plots_to_sort"
    outname="%s/fft_range_%s_n%04d"%(plot_dir,sim,frame)
    bt.plot_brunt(ftool, outname=outname,method='range', fitrange=[4,25])

