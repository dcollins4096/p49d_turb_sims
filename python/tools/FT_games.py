
from GL import *
plt.close('all')

if 0:
    frame=900
    dirname="/data/cb1/Projects/P49_EE_BB/as21_M0.6_MA0.3_64"
    fname1 = "%s/DD%04d/data%04d"%(dirname,frame,frame)
    ds = yt.load(fname1)
    print("Read",ds)
    cg=ds.covering_grid(0,[0.0]*3,[64]*3)

if 0:
    frame=34
    dirname="/data/cb1/Projects/P49_EE_BB/4_1"
    fname1 = "%s/DD%04d/data%04d"%(dirname,frame,frame)
    ds = yt.load(fname1)
    print("Read",ds)
    cg=ds.covering_grid(0,[0.0]*3,[512]*3)

if 0:
    plt.clf()
    plt.imshow( cg['density'].sum(axis=0).v)
    plt.savefig('plots_to_sort/rho1')

import fourier_tools_py3.fourier_filter as Filter
if 0:
    #rho = cg['density'].v

    fft1 = np.fft.fftn( rho )
    ff = Filter.FourierFilter(rho)
    power_1d = np.array([rho[ff.get_shell(bin)].sum() for bin in range(ff.nx)])
    kspace=ff.get_shell_k()

if 0:
    #works
    dirname="/data/cb1/Projects/P49_EE_BB/Products/4_1"
    frame=34
    k3d, rho = dt.dpy('%s/DD%04d.products/power_density.h5'%(dirname,frame), ['k','power'])
    fft1 = np.fft.fftn( rho )
    ff = Filter.FourierFilter(rho)
    power_1d = np.array([rho[ff.get_shell(bin)].sum() for bin in range(ff.nx)])
    kspace=ff.get_shell_k()

    fig,ax=plt.subplots(1,1)
    ax.plot(kspace,power_1d)
    ax.set(xscale='log',yscale='log')
    fig.savefig('plots_to_sort/fft2')

if 0:
    #also works
    frame=34
    dirname="/data/cb1/Projects/P49_EE_BB/4_1"
    fname1 = "%s/DD%04d/data%04d"%(dirname,frame,frame)
    ds = yt.load(fname1)
    cg=ds.covering_grid(0,[0.0]*3,[512]*3)
    rho=cg['density'].v
    fft1 = np.fft.fftn( rho )
    power=fft1*np.conjugate(fft1)
    ff = Filter.FourierFilter(power)
    power_1d = np.array([power[ff.get_shell(bin)].sum() for bin in range(ff.nx)])
    kspace=ff.get_shell_k()

    fig,ax=plt.subplots(1,1)
    ax.plot(kspace,power_1d)
    ax.set(xscale='log',yscale='log')
    fig.savefig('plots_to_sort/fft3')

if 0:
    #small thing for games.
    print('test')
    if 0:
        frame=900
        dirname="/data/cb1/Projects/P49_EE_BB/as21_M0.6_MA0.3_64"
        Nz = 64
    if 1:
        frame=34
        dirname="/data/cb1/Projects/P49_EE_BB/3_2"
        Nz = 512

    fname1 = "%s/DD%04d/data%04d"%(dirname,frame,frame)
    ds = yt.load(fname1)

    cg=ds.covering_grid(0,[0.0]*3,[Nz]*3)
    rho=cg['density'].v
    fft1 = np.fft.fftn( rho )
    power=fft1*np.conjugate(fft1)
    ff = Filter.FourierFilter(power)
    power_1d = np.array([power[ff.get_shell(bin)].sum() for bin in range(ff.nx)])
    Nzones1 = np.array([ff.get_shell(bin).sum() for bin in range(ff.nx)])
    kspace=ff.get_shell_k()

    rho2 = rho.sum(axis=0)
    fft2 = np.fft.fftn( rho2 )
    power2=fft2*np.conjugate(fft2)
    ff2 = Filter.FourierFilter(power2)
    power_1d2 = np.array([power2[ff2.get_shell(bin)].sum() for bin in range(ff2.nx)])
    Nzones2 = np.array([ff2.get_shell(bin).sum() for bin in range(ff2.nx)])
    kspace2=ff2.get_shell_k()

if 0:
    #the good one.
    fig,ax=plt.subplots(1,1)
    ax.set(xscale='log',yscale='log')
    ax.plot(kspace,power_1d/Nzones1,label='3d')
    fig.savefig('plots_to_sort/arg')
    ax.plot(kspace2,power_1d2/Nzones2,label='2d')
    ax.legend(loc=0)
    fig.savefig('plots_to_sort/fft4')


if 1:
    fig,ax=plt.subplots(1,1)
    ax.set(xscale='log',yscale='log')
    Nzones1b=4*np.pi*kspace**2
    ax.plot(kspace,power_1d/Nzones1b,label='3d')
    fig.savefig('plots_to_sort/arg')
    Nzones2b = np.pi*kspace2
    ax.plot(kspace2,power_1d2/Nzones2b,label='2d')
    ax.legend(loc=0)
    fig.savefig('plots_to_sort/fft4b')

if 1:
    ax.clear()
    ax.plot( kspace, Nzones1)
    Kvol=4*np.pi*kspace**2
    fac=Nzones1[2]/Kvol[2]
    print(fac)
    ax.plot( kspace, Kvol*fac,c='r')
    ax.plot(kspace2, Nzones2)
    ax.set(yscale='log')

    fig.savefig('plots_to_sort/Nzones')


if 0:
    #small thing for games.
    print('test')
    frame=900
    dirname="/data/cb1/Projects/P49_EE_BB/as21_M0.6_MA0.3_64"
    fname1 = "%s/DD%04d/data%04d"%(dirname,frame,frame)
    ds = yt.load(fname1)
    cg=ds.covering_grid(0,[0.0]*3,[64]*3)
    rho=cg['density'].v
    rho = rho.sum(axis=0)
    fft1 = np.fft.fftn( rho )
    power=fft1*np.conjugate(fft1)
    ff = Filter.FourierFilter(power)
    power_1d = np.array([power[ff.get_shell(bin)].sum() for bin in range(ff.nx)])
    kspace=ff.get_shell_k()

    #fig,ax=plt.subplots(1,1)
    ax.plot(kspace,power_1d)
    ax.set(xscale='log',yscale='log')
    fig.savefig('plots_to_sort/fft4')

if 0:
    fig,ax=plt.subplots(1,1)
    ax.plot(kspace,power_1d)
    fig.savefig('plots_to_sort/fft1')




