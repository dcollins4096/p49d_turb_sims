
from dtools.starter1 import *


from tifffile import imread
import tiff_poker as TP
plt.close('all')
reload(TP)

V = TP.viewer()
V.guess_scale()
#V.image(vmin=V.minmax[0],vmax=V.minmax[1])
V.image1()
V.image(a=650,b=1800,c=320,d=1000,vmin=V.minmax[0],vmax=V.minmax[1])
