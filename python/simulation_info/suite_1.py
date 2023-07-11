from GL import *
import simulation
reload(simulation)
import sim_colors


for sim in sim_colors.simlist:
    simulation.sim(sim, data_location=dl.sim_dir_base+sim, product_location=dl.product_dir_base+sim, ms=sim_colors.Ms[sim], ma=sim_colors.Ma[sim],
                   color=sim_colors.color[sim],linestyle=sim_colors.linestyle[sim],marker=sim_colors.marker[sim],
                   framelist=sim_colors.framedict[sim])
