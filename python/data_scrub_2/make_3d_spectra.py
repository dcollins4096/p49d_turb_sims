from GL import *
from queb3 import powerlaw_fit as plfit

import simulation

def make_spec(simlist):
    for nsim,sim in enumerate(simlist):
        this_sim=simulation.corral[sim]


        prefix=None
        pack = queb3.simulation_package(  directory=this_sim.data_location,frames=this_sim.framelist,
                                            product_directory=this_sim.product_location, simname=sim)
        for frame in this_sim.framelist:
            pack.make_spectra(frame)
