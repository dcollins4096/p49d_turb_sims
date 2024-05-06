#
# Rough cut the data from the origin TIFF.
#

from dtools.starter1 import *


from tifffile import imread
import tiff_poker as TP
plt.close('all')
reload(TP)

if 'trim' not in dir():
    trim = {}


#not the most elegant of code.
#a,b,c,d are left, right, bottom, and top in pixels.
for DO in [0,1,2,3,4,5]:
    if DO==0:
        which=0
        a=550;b=1800;c=300;d=1000
        na=0;nb=600;nc=0;nd=1000
        section = 'r0_t1'

    if DO==1:
        which=0
        a=100;b=1700;c=1100;d=1900
        na=0;nb=600;nc=0;nd=1000
        section = 'r0_t2'

    if DO==2:
        which=1
        section = 'r60_t1'
        a=550;b=1800;c=300;d=1000
        na=50;nb=600;nc=0;nd=1000

    if DO==3:
        which=1
        section = 'r60_t2'
        a=100;b=1700;c=1100;d=1900
        na=50;nb=600;nc=0;nd=1000


    if DO==4:
        which=2
        section = 'r120_t1'
        a=550;b=1800;c=300;d=1100
        na=50;nb=600;nc=0;nd=1000

    if DO==5:
        which=2
        section = 'r120_t2'
        a=100;b=1700;c=1100;d=2200
        na=50;nb=600;nc=0;nd=1000

    if 1:
        #plot the regions
        V = TP.viewer(which=which)
        V.guess_scale()
        #V.image(vmin=V.minmax[0],vmax=V.minmax[1])
        V.image1('plots_to_sort/%s_full_image'%section)

    if 1:
        #Get the noise level for trimming
        x,y,noise=V.xtract(a=na,b=nb,c=nc,d=nd)
        vmax = noise.max()
        sigma_n = noise.std()
        x,y,noise=V.xtract_and_image(a=na,b=nb,c=nc,d=nd,vmin=0,vmax=vmax, fname='plots_to_sort/%s_noises'%section,zero=True)

    if 1:
        #auto-trim based on the noise level.
        x,y,rough=V.xtract_and_image(a=a,b=b,c=c,d=d,vmin=V.minmax[0],vmax=V.minmax[1], fname='plots_to_sort/%s_rough'%section,zero=False)
        trim[section] = TP.trimmer(rough, sigma_n=2*sigma_n,fname='plots_to_sort/%s_trim'%section, vmin=V.minmax[0],vmax=V.minmax[1])

if 0:
    fptr=h5py.File('p68_laser/ALL_TRIM1.h5','w')
    for section in trim:
        fptr[section]=trim[section]
    fptr.close()
