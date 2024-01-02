from starter1 import *
import yt
import downsample.volavg as volavg
import filament.hessian as hessian
import tools.pcolormesh_helper as pch
import filament.tools as htools
reload(htools)

projax=0
if 'slice_coord' not in dir():
    slice_coord=12

def mslice(arr,axis):
    sl = [slice(None),slice(None),slice(None)]
    sl[axis]=slice_coord
    b= arr[tuple(sl)]
    return b
def asinh(arr,sigma,axis):
    sl = [slice(None),slice(None),slice(None)]
    sl[axis]=slice_coord
    b= np.abs(arr[tuple(sl)])
    #b=np.abs(b)
    const=0.25
    c=const*np.arcsinh(b/(const*sigma))
    c/=c.max()
    return c
def log1(arr,axis):
    sl = [slice(None),slice(None),slice(None)]
    sl[axis]=slice_coord
    b= arr[tuple(sl)]
    c = np.log(np.abs(b))
    c /= c.max()
    return c
def shift(arr,axis):
    sl = [slice(None),slice(None),slice(None)]
    sl[axis]=slice_coord
    c= arr[tuple(sl)]
    c -= c.min()
    c /= c.max()
    return c

def color_packer(a,b,c):
    return np.stack([a,b,c]).transpose()

def color_vec(a,vec):
    return np.stack([a*vec[0],a*vec[1],a*vec[2]]).transpose()

def color_game_1(estuff,outname=None):
    fig,ax=plt.subplots(2,2,figsize=(8,8))

    sigma = estuff.e1.std()
    a1 = asinh(estuff.e0,sigma,0)
    a2 = asinh(estuff.e1,sigma,0)
    a3 = asinh(estuff.e2,sigma,0)
    zero=np.zeros_like(a1)

    ax[0][0].imshow( color_packer(a1,zero,zero))
    ax[0][1].imshow( color_packer(zero,a2,zero))
    ax[1][0].imshow( color_packer(zero,zero,a3))
    ax[1][1].imshow( color_packer(a1,a2,a3))
    fig.tight_layout()
    fig.savefig(outname)
    plt.close(fig)

def color_game_2(estuff,outname=None):
    fig,ax=plt.subplots(1,1,figsize=(8,8))
    sigma = estuff.e1.std()
    a1 = asinh(estuff.e0,sigma,0)
    a2 = asinh(estuff.e1,sigma,0)
    a3 = asinh(estuff.e2,sigma,0)
    zero=np.zeros_like(a1)
    ax.imshow(color_packer(a1,zero,zero))
    fig.savefig('%s_e0'%outname)
    ax.clear()
    ax.imshow(color_packer(a2,zero,zero))
    fig.savefig('%s_e1'%outname)
    ax.clear()
    ax.imshow(color_packer(a3,zero,zero))
    fig.savefig('%s_e2'%outname)
    plt.close(fig)

def color_game_3(estuff,outname=None):
    fig,ax=plt.subplots(1,3,figsize=(8,4))
    sigma = estuff.e1.std()
    a1 = asinh(estuff.e0,sigma,0)
    a2 = asinh(estuff.e1,sigma,0)
    a3 = asinh(estuff.e2,sigma,0)
    v1=color_vec(a1,[1,0.5,0])
    ax[0].imshow(v1)
    v2=color_vec(a3,[0,0.5,1])
    ax[1].imshow(v2)
    ax[2].imshow(v1+v2)
    #ax[1][0].imshow(packer(rho,rho,rho))
    #ax[1][0].imshow(packer(rho,rho,rho))
    #ax[1][1].imshow(packer(a1,a3,rho))
    fig.savefig('%s_game2'%outname)
    plt.close(fig)

def color_game_4(estuff,outname=None):
    fig,ax=plt.subplots(2,2,figsize=(8,8))
    sigma = estuff.e1.std()
    rho = shift(estuff.arr,0)
    a1 = asinh(estuff.e0,sigma,0)
    a2 = asinh(estuff.e1,sigma,0)
    a3 = asinh(estuff.e2,sigma,0)
    zero=np.zeros_like(a1)
    ax[0][0].imshow(color_packer(rho,zero,zero))
    ax[1][0].imshow(color_packer(zero,a1,zero))
    ax[0][1].imshow(color_packer(zero,zero,a3))
    ax[1][1].imshow(color_packer(rho,a1,a3))
    #ax[1].imshow(color_packer(a1,a2,a3))
    fig.savefig('%s_game4'%outname)
    plt.close(fig)


def color_game_5(estuff,outname=None):
    fig,ax=plt.subplots(2,2,figsize=(8,8))
    sigma = estuff.e1.std()
    rho = shift(estuff.arr,0)
    a1 = asinh(estuff.e0,sigma,0)
    a2 = asinh(estuff.e1,sigma,0)
    a3 = asinh(estuff.e2,sigma,0)
    zero=np.zeros_like(a1)
    ax[0][0].imshow(color_packer(rho,rho,rho))
    ax[1][0].imshow(color_packer(a1,a1,zero))
    ax[0][1].imshow(color_packer(zero,a3,a3))
    #ax[1][1].imshow(color_packer(rho,rho,ro))
    if 0:
        b1 = asinh(np.abs(estuff.e1)+np.abs(estuff.e0)+np.abs(estuff.e2),sigma,0)
        mask=(b1>0.4)*(mslice(estuff.e0,0)<0)*(mslice(estuff.e1,0)<0)
    if 1:
        b1 = asinh(np.abs(estuff.e1)+np.abs(estuff.e0)+np.abs(estuff.e2),sigma,0)
        mask=(b1>0.4)*(mslice(estuff.e0,0)<0)*(mslice(estuff.e1,0)<0)

    ok = rho * mask
    img = color_packer(ok,ok,zero)
    img[~mask.transpose()]=[1,1,1]
    ax[1][1].imshow(img)
    if 0:
        hargs={'histtype':'step','bins':np.linspace(-10,10,1000)}
        ax[1][1].hist( mslice(estuff.e0/sigma,projax).flatten(),**hargs)
        ax[1][1].hist( mslice(estuff.e1/sigma,projax).flatten(),**hargs)
        ax[1][1].hist( mslice(estuff.e2/sigma,projax).flatten(),**hargs)
        ax[1][1].set(yscale='log')

    #ax[1][1].imshow(color_packer(rho,a1,a3))
    #ax[1].imshow(color_packer(a1,a2,a3))
    fig.savefig('%s_game5'%outname)
    plt.close(fig)


def color_game_6(estuff,outname=None):
    fig,ax=plt.subplots(2,2,figsize=(8,8))
    sigma = estuff.e1.std()
    rho = shift(estuff.arr,0)
    e1p=np.zeros_like(estuff.e1)
    e1m=np.zeros_like(estuff.e1)
    e1p[estuff.e1>0]=estuff.e1[estuff.e1>0]
    e1m[estuff.e1<0]=estuff.e1[estuff.e1<0]
    #a1 = asinh(estuff.e0,sigma,0)
    #a2 = asinh(estuff.e1,sigma,0)
    #a3 = asinh(estuff.e2,sigma,0)
    a1p = asinh(e1p,sigma,projax)
    a1m = asinh(e1m,sigma,projax)
    zero=np.zeros_like(a1p)

    #ax[0][0].imshow(color_packer(a1p,a1p,zero))
    #ax[0][1].imshow(color_packer(zero,a1m,a1m))
    #ax[1][1].imshow(color_packer(a1p,a1m,zero))
    v1 = color_vec( a1p, [1,0.5,0])
    v2 = color_vec( a1m, [0,0.5,1])
    ax[0][0].imshow(v1)
    ax[0][1].imshow(v2)
    ax[1][0].imshow(v1+v2)


    fig.savefig('%s_game6'%outname)
    plt.close(fig)


