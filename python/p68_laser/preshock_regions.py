
from dtools.starter1 import *
import power_spectrum as ps

import scipy.stats

from tifffile import imread
import tiff_poker as TP
reload(TP)
plt.close('all')

if 'TNA' not in dir():
    fname = 'p68_laser/TRIM_ALIGN.h5'
    fptr=h5py.File(fname,'r')
    TNA = {}
    for field in fptr:
        TNA[field]=fptr[field][()]
    fptr.close()

if 1:
    preshock_region={}
    postshock_region={}
    preshock_abcd={}
    postshock_abcd={}
    preshock_abcd['r60_t1'] = [100,400,100,400]
    postshock_abcd['r60_t1'] = [600,900,100,400]
    preshock_abcd['r60_t2'] = [100,400,100,400]
    postshock_abcd['r60_t2'] = [600,900,100,400]

    preshock_abcd['r120_t1'] = [100,400,100,400]
    postshock_abcd['r120_t1'] = [600,900,100,400]
    preshock_abcd['r120_t2'] = [100,400,100,400]
    postshock_abcd['r120_t2'] = [600,900,100,400]

    preshock_abcd['r0_t1'] = [100,400,100,400]
    postshock_abcd['r0_t1'] = [600,900,100,400]
    preshock_abcd['r0_t2'] = [100,400,100,400]
    postshock_abcd['r0_t2'] = [600,900,100,400]
    for shot in preshock_abcd:
        print(shot)
        arr=TNA[shot]
        V = TP.viewer(arr=arr)
        a,b,c,d=preshock_abcd[shot]
        preshock_region[shot] = V.xtract_and_image(a=a,b=b,c=c,d=d,vmin=None,vmax=None,fname="plots_to_sort/%s_preshock.png"%shot)
    for shot in postshock_abcd:
        print(shot)
        arr=TNA[shot]
        V = TP.viewer(arr=arr)
        a,b,c,d=postshock_abcd[shot]
        postshock_region[shot] = V.xtract_and_image(a=a,b=b,c=c,d=d,vmin=None,vmax=None,fname="plots_to_sort/%s_postshock.png"%shot)

if 1:
    #denoise and make an image
    plot_dir='plots_to_sort'
    import get_foam as gf
    reload(gf)
    denoise_preshock={}

    for shot in preshock_region:
        NOI = gf.noisy(preshock_region[shot])
        NOI.plot_denoise(kmin=0,nmax=9,outname='%s/%s_denoise'%(plot_dir,shot))
        denoise_preshock[shot] = NOI.rhoback


