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

if 0:

#p2 = np.array([[1,0,0],[0,1,0],[0.5]*3,[0.75,0.75,0],[1,1,1]]).transpose()
#p2=np.array([[0.5]*3,[1,1,1]]).transpose()

    a = np.array([0,0,0])/N
    b = np.array([1,1,1])/N
    a.shape=a.size,1
    b.shape=b.size,1
    p1 = np.array([3,3,3])/N
    total = np.zeros([N,N,N])
    for nfil in np.arange(10):
        print("make filament",nfil)
        a = np.random.random(3)
        a.shape=a.size,1
        b = np.random.random(3)
        b.shape=b.size,1
        
        r1=(lineseg_dist( p2,a,b))
        r1.shape=N,N,N
        g1=np.exp(-r1**2/(2*dx))
        total += g1

    r1.shape=N,N,N
    print(r1)
    plt.close('all')
    fig,axes=plt.subplots(2,2)

    if 0:
        axes[0][0].imshow(r1.sum(axis=0),origin='lower',interpolation='nearest')
        axes[0][1].imshow(r1.sum(axis=1),origin='lower',interpolation='nearest')
        axes[1][0].imshow(r1.sum(axis=2),origin='lower',interpolation='nearest')
        fig.savefig('%s/fake'%plotdir)
    if 0:
        axes[0][0].imshow(g1.sum(axis=0),origin='lower',interpolation='nearest')
        axes[0][1].imshow(g1.sum(axis=1),origin='lower',interpolation='nearest')
        axes[1][0].imshow(g1.sum(axis=2),origin='lower',interpolation='nearest')
        fig.savefig('%s/flake'%plotdir)
    if 1:
        axes[0][0].imshow(total.sum(axis=0),origin='lower',interpolation='nearest')
        axes[0][1].imshow(total.sum(axis=1),origin='lower',interpolation='nearest')
        axes[1][0].imshow(total.sum(axis=2),origin='lower',interpolation='nearest')
        fig.savefig('%s/snake'%plotdir)
#print('w',r1)
#print('x',r2)
#print('wtf')
