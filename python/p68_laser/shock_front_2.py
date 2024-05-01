from dtools.starter1 import *

import shot
reload(shot)

def equal_prob(arr,Nbins,ax):
    v = copy.copy(arr)
    Npoints = v.size//Nbins
    v.sort()
    ind = np.arange(0,v.size,Npoints)
    #ind = np.concatenate([ind,v.size-1])
    edges = v[list(ind)]
    print(edges.size)
    cen = 0.5*(edges[1:]+edges[:-1])
    wid = edges[1:]-edges[:-1]
    hist = 1/wid
    ax.bar(cen,hist,width=wid)

def plot_shock_front_2(device, axes):

    axes[0].plot(device.shot1.x, device.shot1.rhobar, c='r')
    axes[0].plot(device.shot2.x, device.shot2.rhobar, c=[0.5]*4)

    axes[0].plot(device.x_2, device.rhobar_2, c='b')
    v = device.vel.in_units('km/s')
    equal_prob( v.v, 16, axes[1])
    axes[1].set(xlabel=v.units)





    

names = ['r120']
if 'devices' not in dir() or True:
    devices={}

for name in names:
    if name not in devices:
        tmp=shot.device(name)
        tmp.bumper([175,300])
        tmp.compute_velocity([200,450])
        print('dx',tmp.shift_x)
        devices[name]=tmp

if 1:
    fig,axes=plt.subplots(3,5,figsize=(15,12))
    ax0=axes[0][0];ax1=axes[0][1];ax2=axes[0][2]; ax3=axes[0][3]; ax4=axes[0][4]
    ax5=axes[1][0];ax6=axes[1][1];ax7=axes[1][2]; ax8=axes[1][3]; ax9=axes[1][4]
    ax10=axes[2][0];ax11=axes[2][1];ax12=axes[2][2]; ax13=axes[2][3]; ax14=axes[2][4]
    for name in names:
        device = devices[name]
        plot_shock_front_2(device, axes[0])
    fig.tight_layout()
    fig.savefig('plots_to_sort/fronts2.pdf')

