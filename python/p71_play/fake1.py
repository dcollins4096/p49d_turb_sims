from starter1 import *

plt.close('all')


import filament.fake_filament as fil
import filament.hessian as hessian
import filament.tools as htools
reload(fil)
reload(htools)

print('good morning')
plotdir = "%s/PigPen"%os.environ['HOME']

rho = fil.make_random_2(64, 128)

if 'esystem' not in dir() or True:
    #the cut_periodic is necessary for sets that aren't actually periodic.
    #They get bad values.
    esystem=htools.eigen_stuff(rho, cut_periodic=True)

if not esystem.done:
    esystem.do()

if 1:
    htools.hist_evals(esystem,"%s/ev2"%plotdir)
if 1:
    htools.proj(esystem,"%s/proj"%plotdir)


