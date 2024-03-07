
from GL import *
import simulation

import fit_many_spectra as fms
reload(fms)
import simulation_info.all_sims as all_sims
plt.close('all')

simlist = all_sims.lists['suite1']
#simlist = ['half_half']

lefts = [3,4,5]
rights = [24,25,26,27]

if 0:
    fms.fitter(simlist, lefts=lefts, rights=rights)
products = ['density','velocity','magnetic',
            'ClTTx','ClTTy','ClTTz',
            'ClEEx','ClEEy','ClEEz',
            'ClBBx','ClBBy','ClBBz']
#simlist=simlist[0:1]
#products=procuts[0:1]
if 1:
    for P in products:
        fms.plot_spectra(simlist,what='hist', field=P)

if 0:
    for P in products:
        fms.plot_spectra(simlist,what='frame', field=P)
if 0:
    for P in products:
        fms.plot_spectra(simlist,what='lefts', field=P)
        fms.plot_spectra(simlist,what='rights', field=P)

