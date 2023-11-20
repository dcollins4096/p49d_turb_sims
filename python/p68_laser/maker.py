from GL import *
import simulation
reload(simulation)
import simulation_info.all_sims as all_sims



import p1_spectra as p1
reload(p1)

sim_list = all_sims.lists['suite1']

if 'ftool' not in dir():
    ftool=None
ftool=p1.brunt_spectra(sim_list)
#p1.drill(['half_1'])
