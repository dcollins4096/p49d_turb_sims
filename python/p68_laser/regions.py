
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

if 'preshock_region' not in dir():
    preshock_region={}
    preshock_abcd={}
    #abcd = left, right, bottom, top
    preshock_abcd['r60_t1'] = [100,400,100,400]
    preshock_abcd['r60_t2'] = [100,400,100,400]

    preshock_abcd['r120_t1'] = [100,400,100,400]
    preshock_abcd['r120_t2'] = [100,400,100,400]

    preshock_abcd['r0_t1'] = [100,400,100,400]
    preshock_abcd['r0_t2'] = [100,400,100,400]
    for shot in preshock_abcd:
        print(shot)
        arr=TNA[shot]
        V = TP.viewer(arr=arr)
        a,b,c,d=preshock_abcd[shot]
        X,Y,Z = V.xtract(a=a,b=b,c=c,d=d)
        preshock_region[shot] = Z



def image_preshock():
    for shot in preshock_abcd:
        print(shot)
        arr=TNA[shot]
        V = TP.viewer(arr=arr)
        a,b,c,d=preshock_abcd[shot]
        fname = 'plots_to_sort/preshock_%s'%shot
        X,Y,Z = V.xtract_and_image(a=a,b=b,c=c,d=d,vmin=None,vmax=None,fname=fname)


if 0:
    #denoise and make an image
    plot_dir='plots_to_sort'
    import get_foam as gf
    reload(gf)
    denoise_preshock={}

    for shot in preshock_region:
        NOI = gf.noisy(preshock_region[shot])
        NOI.plot_denoise(kmin=0,nmax=9,outname='%s/%s_denoise'%(plot_dir,shot))
        denoise_preshock[shot] = NOI.rhoback

                                         
