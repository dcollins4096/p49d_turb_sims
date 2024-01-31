from starter1 import *
import yt
from downsample import volavg
import fourier_tools_py3.fourier_filter as Filter
import brunt_tools as bt
reload(bt)
plt.close('all')
import slabimg
reload(slabimg)

#DEFINE PLOTDIR
plotdir='%s/PigPen/'%(os.environ['HOME'])

sim='6_2'
frame=50

if 'reg' not in dir():
    reg={}

class slabber():
    def __init__(self, sim, frame):

        self.rho_full, self.rho, self.rho_128 = bt.get_cubes(sim,frame, do_rho_4=True)
        self.sim=sim
        self.frame=frame

        self.ftool = bt.fft_tool(self.rho_128)
        self.ftool.do3()
        self.ftool.do2(projax=0)
    def do_slab(self,projax=0):
        slabimg.slab(self.ftool,projax=0)

if '6_2' not in reg:
    sl = slabber('6_2',50)
    sl.do_slab(projax=0)
    reg['6_2']=sl
if 'half_half' not in reg:
    sl = slabber('half_half',30)
    sl.do_slab(projax=0)
    reg['half_half']=sl

if 0:
    fig,ax=plt.subplots(1,1)
    for sim in ['6_2','half_half']:
        print('Sim',sim)
        outname = "%s/partial_%s"%(plotdir,sim)
        slabimg.number_checker(reg[sim].ftool, axin=ax,label=sim)
        ax.legend(loc=0)
    fig.savefig(outname)

if 1:
    for sim in ['6_2','half_half']:
        outname = "%s/proj_%s"%(plotdir,sim)

        slabimg.plot_fft(reg[sim].ftool,outname=outname)


if 0:
    for sim in ['6_2','half_half']:
        outname = "%s/dissassemble_%s"%(plotdir,sim)
        slabimg.plot(reg[sim].ftool,outname=outname)
    
