

from GL import *

def fig_squares(nx,ny,cbar=True,figwidth=4):


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
    print("LEFT",lefts)
    print("BOTTOMS",bottoms)



    print(px,py)
    print('width %f height %f'%(px,py))

    fig=plt.figure(figsize=(figwidth,figheight))
    print('figheight',figheight,figwidth)
    axes=[]
    for iy in range(ny):
        axes.append([])
        for ix in range(nx):
            print('=====',ix,iy)
            rect = [lefts[ix], bottoms[iy], px, py]
            #rect = [0.1,0.1,0.8,0.8]
            axes[iy].append( plt.axes( rect))
            print(rect)
    cbar_rect = [lefts[-1]+px,Bmargin,cbar_width,1-Tmargin-Bmargin]
    cbar = plt.axes(cbar_rect)
    axes=nar(axes)


    return fig, axes, cbar



