from GL import *
import simulation
import spectra_tools as st

def make_spec(simlist):
    for sim_name in simlist:
        this_sim = simulation.corral[sim_name]
        print(this_sim.all_frames)
        if this_sim.B_nom > 0:
            do_magnetic=True
        else:
            do_magnetic=False
        for frame in this_sim.ann_frames:
            oober = st.short_oober(this_sim.data_location, frame=frame, product_directory=this_sim.product_location, simname=this_sim.name, code=this_sim.code)
            st.MakeDensitySpectra(oober,frame)
            st.MakeVelocitySpectra(oober,frame)
            if do_magnetic:
                st.MakeMagneticSpectra(oober,frame)

def old_make_spec(simlist):
    for nsim,sim in enumerate(simlist):
        this_sim=simulation.corral[sim]


        prefix=None
        pack = queb3.simulation_package(  directory=this_sim.data_location,frames=this_sim.all_frames,
                                            product_directory=this_sim.product_location, simname=sim)
        for frame in this_sim.all_frames:
            print("Spectra on ",sim,frame)
            pack.make_spectra(frame)
