


#
# Make hessian on downsampled cube.
#


from starter1 import *
import yt
import downsample.volavg as volavg
import filament.hessian as hessian
import tools.pcolormesh_helper as pch
import filament.tools as htools
import filament.color_games as cgames
reload(cgames)
reload(htools)

plotdir = "%s/PigPen"%(os.environ['HOME'])


sim = '6_half'
sim = '1_1'
frame = 31



if 'esystem' not in dir():
    print('get')
    rho_full, rho256 = htools.get_cubes(sim,frame)
    print('esystem')

    esystem = htools.eigen_stuff(np.log(rho256))
    esystem.do()
    esystem.doproj(axis=0)

if 1:
    htools.pproj(esystem,"%s/2d_vs_3d_%s_n%04d"%(plotdir,sim,frame))
if 0:
    htools.pproj_cgames(esystem,"%s/derp"%plotdir)
