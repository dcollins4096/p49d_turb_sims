

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
frame = 36



if 'esystem' not in dir():
    print('get')
    rho_full, rho256 = htools.get_cubes(sim,frame)
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

if 0:
    #removed
    htools.color_game(esystem,"%s/turb_color_game"%plotdir)

if 0:
    #3 eigenvalues in a channel map
    cgames.color_game_1(esystem,"%s/slice_game_e0_e1_e2_%s_n%04d"%(plotdir,sim,frame))
if 0:
    #individual eigen vectors
    cgames.color_game_2(esystem,"%s/slice_game_%s_n%04d"%(plotdir,sim,frame))

if 0:
    #just e0 and e2
    cgames.color_game_3(esystem,"%s/slice_game_e0_e2_%s_n%04d"%(plotdir,sim,frame))

if 0:
    cgames.color_game_4(esystem,"%s/slice_game_density_%s_n%04d"%(plotdir,sim,frame))

if 1:
    cgames.color_game_5(esystem,"%s/slice_game_cutoff_%s_n%04d"%(plotdir,sim,frame))
if 1:
    cgames.color_game_6(esystem,"%s/slice_game_e2pos_%s_n%04d"%(plotdir,sim,frame))





