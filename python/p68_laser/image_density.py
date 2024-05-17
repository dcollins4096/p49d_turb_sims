
from dtools.starter1 import *


from tifffile import imread
import tiff_poker as TP
reload(TP)
import p68_laser.numbers as numbers
reload(numbers)

if 'TNA' not in dir():
    fname = 'p68_laser/TRIM_ALIGN.h5'
    fptr=h5py.File(fname,'r')
    TNA = {}
    for field in fptr:
        TNA[field]=fptr[field][()]
    fptr.close()



fig, axes = plt.subplots(3,2,figsize=(12,12))
for nr,pair in enumerate([['r0_t1','r0_t2'],['r60_t1','r60_t2'],['r120_t1','r120_t2']]):
    for nc in [0,1]:
        q = pair[nr]


        axes[nr,nc].imshow(


