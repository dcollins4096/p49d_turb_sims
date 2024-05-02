from GL import *

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

if 1:
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
