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
    ax1.hist(vel.v)
    ax1.set(xlabel=vel.units)
    return dx




def plot_front(name, ax,ax2=None, **kwargs):
    #Q = phys.useful_values_take1(name)
    Q = phys.useful_values_take2(name)
    print('fart')
    #print("Negative values:", (Q<=0).sum())

    sl = slice(200,450)
    std = Q[sl,:].std(axis=0)
    Qbar = Q[sl,:].mean(axis=0)
    ax.plot(Q[sl,:].transpose(), linewidth=0.2,c=[0.5,0.5,0.5,0.1])
    ax.plot(Qbar, **kwargs)
    if ax2:
        xax = np.arange(0,Q.shape[0])
        ax2.imshow(Q)
        for line in xax[sl]:
            ax2.axhline(line,c=[0.5,0.5,0.5,0.1])



if 1:
    fig,axes=plt.subplots(3,5,figsize=(15,12))
    ax0=axes[0][0];ax1=axes[0][1];ax2=axes[0][2]; ax3=axes[0][3]; ax4=axes[0][4]
    ax5=axes[1][0];ax6=axes[1][1];ax7=axes[1][2]; ax8=axes[1][3]; ax9=axes[1][4]
    ax10=axes[2][0];ax11=axes[2][1];ax12=axes[2][2]; ax13=axes[2][3]; ax14=axes[2][4]
    plot_front('r0_t1'  ,ax0 ,ax1 ,c='r')
    plot_front('r0_t2'  ,ax0 ,ax2 ,c='g')
    plot_front('r60_t1' ,ax5 ,ax6 ,c='r')
    plot_front('r60_t2' ,ax5 ,ax7 ,c='g')
    plot_front('r120_t1',ax10,ax11,c='r')
    plot_front('r120_t2',ax10,ax12,c='g')

    print('compute velocities (takes a sec)')
    dx0=  diff('r0',   [410,600],ax3, ax4  )
    dx60= diff('r60',  [400,600],ax8, ax9  )
    dx120=diff('r120', [400,600],ax13, ax14  )

    fig.tight_layout()
    fig.savefig('plots_to_sort/shock_fronts_b')

