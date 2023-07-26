

from GL import *

def fig_squares(nx,ny,cbar=True,figwidth=4, verbose=False):


    Lmargin=0.03
    Rmargin=0.15
    Tmargin=0.03
    Bmargin=0.03
    cbar_width=0.03
    px = (1.-Lmargin-Rmargin-cbar_width)/nx
    py = (1.-Bmargin-Tmargin)/ny

    lefts = [Lmargin + px*n for n in range(nx)]
    bottoms = [Bmargin + py*n for n in range(ny)]

    figheight = figwidth*px/py
    if verbose:
        print("LEFT",lefts)
        print("BOTTOMS",bottoms)
        print(px,py)
        print('width %f height %f'%(px,py))




    fig=plt.figure(figsize=(figwidth,figheight))
    if verbose:
        print('figheight',figheight,figwidth)
    axes=[]
    for iy in range(ny):
        axes.append([])
        for ix in range(nx):
            rect = [lefts[ix], bottoms[iy], px, py]
            axes[iy].append( plt.axes( rect))
    cbar_rect = [lefts[-1]+px,Bmargin,cbar_width,1-Tmargin-Bmargin]
    cbar = plt.axes(cbar_rect)
    axes=nar(axes)[::-1,:]


    return fig, axes, cbar



