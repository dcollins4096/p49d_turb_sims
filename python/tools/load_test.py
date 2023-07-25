


from GL import *

import simulation
import simulation_info.all_sims

def loader(simlist):
    plt.close('all')
    for this_simname in simlist:
        print('slopes',this_simname)
        this_sim = simulation.corral[this_simname]
        #this_sim.read_all_spectra()
        #this_sim.fit_all_spectra()
        this_sim.load()
        #print(this_sim.quan_time['time'].size)
        #this_sim.read_pdfs(['magnetic_field_strength'])
        #print(this_sim.pdfs)



loader(['4_1'])
