from starter1 import *


def lineseg_dist(p, a, b):

    # normalized tangent vector
    d = np.divide(b - a, np.linalg.norm(b - a))

    # signed parallel distance components
    s = np.dot((a - p).transpose(), d)
    #t = np.dot((p - b).transpose(), d)
    # clamped parallel distance
    #h = np.maximum.reduce([s, t, np.zeros_like(t)])
    # perpendicular distance component
    #print(d.shape)
    c = np.cross(p - a, d, axisa=0,axisb=0)

    r1= np.linalg.norm(c,axis=1)
    #r2= np.hypot(h, np.linalg.norm(c))
    return r1

def make_random_1(Nfil=1,Nside=16):
    dx = 1/Nside
    x,y,z=np.mgrid[0:1:dx,0:1:dx,0:1:dx]
    x=x.flatten()
    y=y.flatten()
    z=z.flatten()
    p2 = np.stack([x,y,z])
    p2.transpose()
    total = np.zeros([Nside,Nside,Nside])
    for nfil in np.arange(Nfil):
        print("make filament",nfil)
        a = np.random.random(3)
        a.shape=a.size,1
        b = np.random.random(3)
        b.shape=b.size,1
        
        r1=(lineseg_dist( p2,a,b))
        r1.shape=Nside,Nside,Nside
        g1=np.exp(-r1**2/(2*dx))
        total += g1
    return total

def make_random_2(Nfil=1,Nside=16,sigma=None):
    dx = 1/Nside
    x,y,z=np.mgrid[0:1:dx,0:1:dx,0:1:dx]
    x=x.flatten()
    y=y.flatten()
    z=z.flatten()
    p2 = np.stack([x,y,z])
    p2.transpose()
    total = np.zeros([Nside,Nside,Nside])
    if sigma==None:
        sigma=0.1*dx
    for nfil in np.arange(Nfil):
        print("make filament",nfil)
        a = np.random.random(3)
        a.shape=a.size,1
        b = np.random.random(3)
        b.shape=b.size,1
        
        r1=(lineseg_dist( p2,a,b))
        r1.shape=Nside,Nside,Nside
        g1=np.exp(-r1**2/(2*sigma))
        total += g1
    return total

