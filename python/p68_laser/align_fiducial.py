

from dtools.starter1 import *


from tifffile import imread
import tiff_poker as TP
plt.close('all')
reload(TP)
import aligner
reload(aligner)
plt.close('all')
if 'trimmed' not in dir():
    fptr=h5py.File('p68_laser/ALL_TRIM1.h5','r')
    trimmed={}
    for section in fptr:
        trimmed[section]=fptr[section][()]
    fptr.close()

#aligner.align(trimmed['r120_t1'], trimmed['r120_t2'], name_base='r120', xrange=[290,350])
aligner.align(trimmed['r60_t1'], trimmed['r60_t2'], name_base='r60', xrange=[300,340], x_override=0)
#aligner.align(trimmed['r0_t1'], trimmed['r0_t2'], name_base='r0', xrange=[300,330])
