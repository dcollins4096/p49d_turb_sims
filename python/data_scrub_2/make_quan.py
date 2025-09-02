from GL import *

import simulation as simulation
reload(simulation)
import compute_avg_quantities as comp_avg
reload(comp_avg)

def comp_all(simlist):
    for sim_name in simlist:
        this_sim = simulation.corral[sim_name]
        print(this_sim.all_frames)
        for frame in this_sim.ann_frames:
            if this_sim.code == 'Enzo':
                comp_avg.make_quan(this_sim.data_location,frame,out_directory=this_sim.product_location,sim=sim_name, clobber=False, do_magnetic=this_sim.do_magnetic )
            elif this_sim.code == 'Athena':
                comp_avg.make_quan_athena(this_sim.data_location,frame,out_directory=this_sim.product_location,sim=sim_name, clobber=False, do_magnetic=this_sim.do_magnetic )


def comp_bulk(simlist):
    for sim_name in simlist:
        this_sim = simulation.corral[sim_name]
        for frame in this_sim.all_frames:
            comp_avg.bulk_viscosity_estimate(this_sim.data_location,frame,out_directory=this_sim.product_location,sim=sim_name, clobber=False )
def energy_cleaner(simlist):
    for sim_name in simlist:
        this_sim = simulation.corral[sim_name]
        for frame in this_sim.all_frames:
            outname = "%s/DD%04d.products/data%04d.AverageQuantities.h5"%(this_sim.product_location,frame,frame)
            fptr = h5py.File(outname,'r+')
            if 'Ekin' in fptr:
                print('yes', frame)
                del fptr['Ekin']
            else:
                print('No Ekin')
            fptr.close()


def comp_Edot(simlist):
    total = 0
    for sim_name in simlist:
        this_sim = simulation.corral[sim_name]
        for frame in this_sim.all_frames:
            total += 1

    done = 0
    import time
    start_time = time.time()

    for sim_name in simlist:
        this_sim = simulation.corral[sim_name]
        for frame in this_sim.all_frames:
            comp_avg.make_edot_faster(this_sim.data_location,frame,out_directory=this_sim.product_location,sim=sim_name, clobber=False )
            tnow = time.time()
            done += 1
            dt = tnow-start_time
            avg_rate = dt/done
            time_left = (total-done)*avg_rate/60
            print( "Finished %d/%d, %f seconds ellapsed = %f minutes Remaining %f"%(done, total, dt, dt/60, time_left))

def comp_Ekin(simlist):
    total = 0
    for sim_name in simlist:
        this_sim = simulation.corral[sim_name]
        for frame in this_sim.all_frames:
            total += 1

    done = 0
    import time
    start_time = time.time()

    for sim_name in simlist:
        this_sim = simulation.corral[sim_name]
        for frame in this_sim.all_frames:
            comp_avg.make_ekin(this_sim.data_location,frame,out_directory=this_sim.product_location,sim=sim_name, clobber=False )
            tnow = time.time()
            done += 1
            dt = tnow-start_time
            print( "Finished %d/%d, %f seconds ellapsed = %f minutes"%(done, total, dt, dt/60))
