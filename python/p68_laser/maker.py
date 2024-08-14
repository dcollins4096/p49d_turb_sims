from GL import *

import simulation
reload(simulation)
import simulation_info.all_sims as all_sims
#import old_brunt_tools as bt
import dtools.math.brunt_tools as bt
reload(bt)


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

if 1:
    sim_list = all_sims.lists['suite1']

    import all_brunt
    reload(all_brunt)
    all_brunt.plot_all_brunt(sim_list,projax=1)
    all_brunt.plot_sigmas(sim_list,projax=1)

if 0:
    N = 128
    alphaT=-1.5
    alphaV=alphaT-2
    kmin=2
    kmax=-2
    Q = bt.fake_powerlaw(N,alphaT,kmin,kmax)
    ftool=bt.fft_tool(Q)
    ftool.do3()
    ftool.do2(projax=0)
    bt.plot_brunt(ftool,method='full',outname='%s/fake_powerlaw'%plotdir)

    plt.close('all')
    fig,axes=plt.subplots(1,2)
    ax0=axes[0];ax1=axes[1]
    M1 = ftool.ps2.power>1e-16
    M2 = M1
    M3 = ftool.ps3.power>1e-16
    kp2d=(ftool.ps2.kcen*ftool.ps2.power)
    ax0.plot( ftool.ps2.kcen[M1], ftool.ps2.power[M1], c='r')
    ax0.plot( ftool.ps3.kcen[M3], ftool.ps3.power[M1], c='g')
    ax0.plot( ftool.ps2.kcen[M1], 2*kp2d[M1],c='b')
    ax0.set(xscale='log',yscale='log')
    import dtools.math.equal_probability_binner as epb
    RRR=kp2d/ftool.ps3.power
    hist, cen, wid=epb.equal_prob( RRR[M1], 16, ax=ax1)
    print( cen[ np.argmax(hist)])
    #TheX, TheY = ftool.ps2.kcen[M1], RRR[M1]
    #pfit = np.polyfit( np.log10(TheX), np.log10(TheY),1)
    #print(pfit)
    #ax1.plot( TheX, TheY)
    #ax1.set(xscale='log', yscale='log')
    fig.savefig('%s/play'%plotdir)

