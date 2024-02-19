

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

TA_120_1, TA_120_2=aligner.align(trimmed['r120_t1'], trimmed['r120_t2'], name_base='r120', xrange=[290,350])
TA_60_1, TA_60_2=aligner.align(trimmed['r60_t1'], trimmed['r60_t2'], name_base='r60', xrange=[300,340], x_override=0)
TA_0_1, TA_0_2=aligner.align(trimmed['r0_t1'], trimmed['r0_t2'], name_base='r0', xrange=[300,330])

if 1:
    trimname = 'p68_laser/TRIM_ALIGN.h5'
    if os.path.exists(trimname):
        print("File exists, skipping", trimname)
    else:
        fptr=h5py.File(trimname,'w')
        fptr['r120_t1'] = TA_120_1
        fptr['r120_t2'] = TA_120_2
        fptr['r60_t1'] = TA_60_1
        fptr['r60_t2'] = TA_60_2
        fptr['r0_t1'] = TA_0_1
        fptr['r0_t2'] = TA_0_2
        fptr.close()
