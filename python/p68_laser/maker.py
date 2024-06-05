from GL import *

import simulation
reload(simulation)
import simulation_info.all_sims as all_sims


if 0:
    import simulation
    reload(simulation)
    import simulation_info.all_sims as all_sims

    import p1_spectra as p1
    reload(p1)

    sim_list = all_sims.lists['suite1']

    if 'ftool' not in dir():
        ftool=None
    sim_list=['6_2']
    ftool=p1.brunt_spectra(sim_list)
    #p1.drill(['half_1'])

if 0:
    import regions
    reload(regions)
    shot=['r60_t1']#, 'r60_t1','r120_t1']
    regions.image_zero(shot)
    print(regions.get_zero(shot[0]))
    #zero = regions.zero_region[shot[0]]

if 0:
    import horizontal_distance as horiz
    reload(horiz)
    horiz.try1(method=1,fname='test1')
    stuff=horiz.try1(method=2,fname='test2')
    #print(stuff)

if 0:
    import data_scrub_2.make_small_rho as msr
    reload(msr)
    msr.make(sim_list)

if 0:
    sim_list = all_sims.lists['suite1']
    import all_brunt
    reload(all_brunt)
    #all_brunt.plot(sim_list)
    all_brunt.plot_sigmas(sim_list)

if 0:
    N = 128
    alphaT=-1.5
    alphaV=alphaT-2
    kmin=2
    kmax=-2
    Q = bt.fake_powerlaw(N,alphaT,kmin,kmax)

