
from dtools.starter1 import *

from tifffile import imread
import tiff_poker as TP


if 0:
    #Coarse Cut from raw radiograph
    #the import will extract the two frames from the raw radiograph.
    #Will save to h5 if asked.
    import regions_coarse

if 0:
    #Trim off noise at the edges.
    #Align to the point of the 2nd fiducial.
    import regions_align


