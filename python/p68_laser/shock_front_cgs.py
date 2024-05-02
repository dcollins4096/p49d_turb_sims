from dtools.starter1 import *

import regions
reload(regions)
import physical_values as phys
reload(phys)

import horizontal_distance as horz
reload(horz)
#inter=horz.try1()

def diff(base,rng,ax0,ax1):
    shot1 = '%s_t1'%base
    shot2 = '%s_t2'%base
    F1 = phys.useful_values_take2(shot1)
    F2 = phys.useful_values_take2(shot2)
    X1 = phys.get_x(shot1)
    X2 = phys.get_x(shot2)
    sl = slice(200,450)

    front1 = F1[sl,:].mean(axis=0)
    front2 = F2[sl,:].mean(axis=0)

    ok = slice(rng[0],rng[1])
    S1 = front1[ok]
    S2 = front2[ok]
    intersections = horz.ho(S2,S1)
    I = nar(intersections)
    dx = I[:,1]-I[:,0]

    #S1 = [4, 1, 2, 7, 8, 8, 6, 11, 7, 10, 11, 15, 14, 14, 13, 17, 17, 21, 22, 20]
    #S2 = [3, 0, 1, 6, 3, 6, 9, 11, 8, 8, 11, 15, 14, 15, 17, 14, 18, 17, 18, 20]
    #print(S1,S2)
    #horz.try2(S2,S1)

    ax0.plot(X1[ok],front1[ok])
    ax0.plot(X2[ok],front2[ok])
    vel = phys.pixel_to_velocity(dx)
    vel.sort()
    N = 16
    ax1.hist(vel.v)
    ax1.set(xlabel=vel.units)
    return vel


import shot
reload(shot)

def plot_front(name, ax,ax2=None, **kwargs):
    #Q = phys.useful_values_take1(name)
    Q = phys.useful_values_take2(name)
    #print("Negative values:", (Q<=0).sum())

    sl = slice(200,400)
    std = Q[sl,:].std(axis=0)
    Qbar = Q[sl,:].mean(axis=0)
    ax.plot(Q[sl,:].transpose(), linewidth=0.2,c=[0.5,0.5,0.5,0.1])
    ax.plot(Qbar, **kwargs)
    if ax2:
        xax = np.arange(0,Q.shape[0])
        ax2.imshow(Q)
        for line in xax[sl]:
            ax2.axhline(line,c=[0.5,0.5,0.5,0.1])
    return Q[sl,:]



if 'vel60' not in dir():
    fig,axes=plt.subplots(3,5,figsize=(15,12))
    ax0=axes[0][0];ax1=axes[0][1];ax2=axes[0][2]; ax3=axes[0][3]; ax4=axes[0][4]
    ax5=axes[1][0];ax6=axes[1][1];ax7=axes[1][2]; ax8=axes[1][3]; ax9=axes[1][4]
    ax10=axes[2][0];ax11=axes[2][1];ax12=axes[2][2]; ax13=axes[2][3]; ax14=axes[2][4]
    rho_0_1=plot_front('r0_t1'  ,ax0 ,ax1 ,c='r')
    rho_0_2=plot_front('r0_t2'  ,ax0 ,ax2 ,c='g')
    rho_60_1=plot_front('r60_t1' ,ax5 ,ax6 ,c='r')
    rho_60_2=plot_front('r60_t2' ,ax5 ,ax7 ,c='g')
    rho_120_1=plot_front('r120_t1',ax10,ax11,c='r')
    rho_120_2=plot_front('r120_t2',ax10,ax12,c='g')

    print('compute velocities (takes a sec)')
    vel0=  diff('r0',   [410,600],ax3, ax4  )
    vel60= diff('r60',  [400,600],ax8, ax9  )
    vel120=diff('r120', [400,600],ax13, ax14  )

    fig.tight_layout()
    fig.savefig('plots_to_sort/shock_fronts_b')
if 1:
    fig,axes=plt.subplots(3,3, figsize=(12,12))
    vels = [vel0,vel60,vel120]
    rho_1 = [rho_0_1,rho_60_1,rho_120_1]
    rho_2 = [rho_0_2,rho_60_2,rho_120_2]


    Nbins=16
    for n,v in enumerate(vels):
        Npoints = v.size//Nbins

        v = copy.copy(v)
        v.sort()
        ind = np.arange(0,v.size,Npoints)
        #ind = np.concatenate([ind,v.size-1])
        edges = v[list(ind)]
        print(edges.size)
        cen = 0.5*(edges[1:]+edges[:-1])
        wid = edges[1:]-edges[:-1]
        hist = 1/wid
        ax=axes[0][n]
        ax.bar(cen,hist,width=wid)
        max_bin = np.argmax(hist)
        sigma_v = cen[max_bin]
        ax.text(0.5,0.5,r'$\sigma_v=%0.2f\rm{km/s}$'%(sigma_v), transform=ax.transAxes)


        ax=axes[1][n]
        rho_1_bar=rho_1[n].mean(axis=0)
        rho_2_bar=rho_2[n].mean(axis=0)
        ax.plot(rho_1_bar)
        ax.plot(rho_2_bar)

        ax2=axes[2][n]
        rng = {0:[175,300],1:[175,300],2:[175,300]}[n]

        ax.axvline(rng[0])
        ax.axvline(rng[1])

        ax2.plot(rho_1_bar[slice(*rng)])
        ax2.plot(rho_2_bar[slice(*rng)])




        
    fig.tight_layout()
    fig.savefig('plots_to_sort/EP.pdf')



