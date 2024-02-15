from tifffile import imread
import matplotlib as mpl
import matplotlib.pyplot as plt
import copy
import numpy as np
plt.close('all')

all_data = imread('TD_TC090-124_HGXD_IMAGE_N220712-002-999_DROOP_CORR_422421128478532_20220907115837893.tif')

if 1:
    #get extents
    s1 = all_data[250:1100,500:]
    vmin = s1.min()
    vmin=0
    vmax = s1.max()
    norm1 = mpl.colors.Normalize(vmin=vmin,vmax=vmax)

    XX = np.arange( all_data.shape[0])
    YY = np.arange( all_data.shape[1])
    X1, Y1 = np.meshgrid(XX,YY)

if 0:
    fig,ax=plt.subplots(1,1)
    for ypos in [0.1,0.2,0.3,0.4,0.5]:
        Y = s1[int(ypos*s1.shape[0]),:]
        ax.plot(Y)
    fig.savefig('1d')


if 1:
    #get extents
    fig,ax=plt.subplots(1,1)
    ax.imshow(s1[500:700,400:600],origin='lower')
    fig.savefig('fiducial')

if 1:
    #get extents
    fig,ax=plt.subplots(1,1)
    ax.imshow(s1[300:500,400:600],origin='lower')
    fig.savefig('noise')

if 1:
    #get extents
    fig,ax=plt.subplots(1,1)
    s3=s1[200:600,1000:1200]
    #ax.imshow(s3,origin='lower',norm=norm)
    ax.imshow(s3,origin='lower')
    fig.savefig('signal1')


if 0:
    #map1.
    fig,ax=plt.subplots(1,1)
    norm=mpl.colors.Normalize(vmin=0, vmax=vmax)
    #ax.imshow(all_data,norm=norm)
    ax.pcolormesh(X1,Y1,all_data,norm=norm)
    fig.savefig('all_0')

if 0:
    fig,ax=plt.subplots(1,1)
    ok = (X1 < 2000)*(X1<1000)*(Y1<1500)
    ax.pcolormesh(X1[ok],Y1[ok],all_data[ok],norm=norm)
    #fig.savefig('shot_2')


if 0:
    s1 = all_data[250:1100,500:]
    vmin = s1.min()
    vmin=0
    vmax = s1.max()
    mpl.image.imsave('full_uniform.png', all_data, cmap='viridis', vmin=vmin,vmax=vmax)
    mpl.image.imsave('f1.png',s1,cmap='viridis', vmin=0)

    plt.clf()
    plt.hist(s1.flatten(),bins=100)
    plt.savefig('s1_flatten')


if 0:
    fig,axes=plt.subplots(1,2)
    ax0=axes[0];ax1=axes[1]

    ax0.hist(numpy_data.flatten())
    norm = mpl.colors.Normalize(vmin=0,vmax=numpy_data.max()/4)
    cmap_name='viridis'
    cmap = copy.copy(mpl.cm.get_cmap(cmap_name))
    cmap.set_under('w')
    ax1.imshow(numpy_data,norm=norm,cmap=cmap)
    fig.savefig('thing1')

