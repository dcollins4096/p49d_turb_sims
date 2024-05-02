from dtools.starter1 import *

import shot
reload(shot)
import equal_probability_binner as ep


def plot_shock_front_2(device, axes):

    axes[0].plot(device.shot1.x, device.shot1.rhobar, c='r')
    axes[0].plot(device.shot2.x, device.shot2.rhobar, c=[0.5]*4)

    axes[0].plot(device.x_2, device.rhobar_2, c='b')
    v = device.vel_dist.in_units('km/s')
    hist,cent=ep.equal_prob( v.v, 16, axes[1])
    axes[1].set(xlabel=v.units)
    most_prob = np.argmax(hist)
    vshock = cent[most_prob]/device.cs

    #sigma_v = 0.5*(60/45)**0.33*device.atwood*vshock/device.cs

    axes[1].text(0.5,0.75,r'$v_s=%0.2f$'%(vshock), transform=axes[1].transAxes)

names = ['r60']#,'r120', 'r0']
if 'devices' not in dir() or True:
    devices={}
if 'r0' in  devices:
    del devices['r0']
vel_cut = {'r120':[325,450], 'r60':[325,450], 'r0':[325,450]}
means = {'r60':[325,375, 650]}
for name in names:
    if name not in devices:
        tmp=shot.device(name)
        #tmp.image(fname = 'image_shot_%s'%name)
        x_off=None
        if name == 'r0':
            shift_60 = devices['r60'].shift_x
            shift_120 = devices['r120'].shift_x
            x_off = 0.5*(shift_60+shift_120)
        tmp.bumper([175,300],fix_shift_x=x_off,fname='bumper_%s.pdf'%name, do_plot=True )
        tmp.compute_velocity( vel_cut[name],fname = 'get_horiz_%s'%name)
        tmp.atwood( mean_density=means[name], fname = 'atwood_%s'%name)
        print('dx',tmp.shift_x)
        devices[name]=tmp

if 1:
    fig,axes=plt.subplots(3,5,figsize=(15,12))
    ax0=axes[0][0];ax1=axes[0][1];ax2=axes[0][2]; ax3=axes[0][3]; ax4=axes[0][4]
    ax5=axes[1][0];ax6=axes[1][1];ax7=axes[1][2]; ax8=axes[1][3]; ax9=axes[1][4]
    ax10=axes[2][0];ax11=axes[2][1];ax12=axes[2][2]; ax13=axes[2][3]; ax14=axes[2][4]
    #plot_shock_front_2(devices['r0'], axes[0])
    plot_shock_front_2(devices['r60'], axes[1])
    #plot_shock_front_2(devices['r120'], axes[2])
    fig.tight_layout()
    fig.savefig('plots_to_sort/fronts2.pdf')

