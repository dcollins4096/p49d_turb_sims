"""
Make many queb data products.

*  queb3.simulation_package points to a simulation.  
   BoxSize is in units of 64 zones.
*  
      
"""
from GL import *
import queb3
reload(queb3)
import sim_colors
plot_dir =  "/home/dccollins/PigPen"

#a thing that describes the simulation
for sim in sim_colors.simlist:
    frames=sim_colors.framedict[sim]
    sim_dir = sim_colors.cloudbreak_base + "/" + sim
    product_dir = sim_colors.cloudbreak_base + "/Products/" + sim
    pack = queb3.simulation_package( directory=sim_dir,frames=frames,prefix=sim, product_directory=product_dir)
#produce all QUEB products.
    #pack.EBall()
    for frame in frames:
        #makes 3d spectra
        pack.make_spectra(frame=frame)
        #read Q&U.  Compute E&B.
        this_proj=pack.read_queb(frame=frame,ax='x')
        this_proj.compute_harmonic_products()
        #read 3d spectra    
        this_proj.read_spectra(frame)
        #here we will make improvements.
        #Also here we can vary the fit range.
        #Further, we can average this_proj.ClEE over several frames to 
        #get a more meaningful result.
        fitrange=this_proj.determine_fit_range()  #something better 
        #do the fit.  Not much fancy here
        slopes=this_proj.fit_eb_slopes(fitrange=fitrange)
        #this plotting package can be generalized.
        this_proj.prefix=sim
        pack.plot_eb(this_proj, fname='%s/%s_x_n%04d.png'%(plot_dir,sim,frame),slopes=slopes)
        pack.plot_many_spectra(this_proj)
