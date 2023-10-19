

#
# Make hessian on downsampled cube.
#


from starter1 import *
import downsample.volavg as volavg
import filament.hessian as hessian
import tools.pcolormesh_helper as pch
import filament.tools as htools
reload(htools)
import yt

plotdir = "%s/PigPen"%(os.environ['HOME'])


sim = '4_1'
frame = 35



if 'ds' not in dir():
    print('load and cg')
    ds = yt.load("/data/cb1/Projects/P49_EE_BB/%s/DD%04d/data%04d"%(sim, frame, frame))
    print('get cg')
    cg = ds.covering_grid(0, [0.0]*3, [512]*3)
    dds = cg.dds
    rho_full = cg["density"].v
    print('volavg')
    rho = volavg.volavg(rho_full,rank=3,refine_by=2)
    rho256=rho
    print('esystem')
    esystem = htools.eigen_stuff(np.log(rho256))
    esystem.do()


if 0:
    htools.hist_evals(esystem, "%s/turb_evals"%plotdir)

if 0:
    htools.proj(esystem,"%s/turb_proj"%plotdir)

if 0:
    htools.doslice(esystem,"%s/turb_slice"%plotdir)

if 0:
    htools.contour(esystem,"%s/turb_contour"%plotdir)

if 0:
    htools.correlations(esystem,"%s/turb_cor"%plotdir)
if 1:
    htools.color_game(esystem,"%s/turb_color_game"%plotdir)



