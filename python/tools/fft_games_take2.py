from starter1 import *
plt.close('all')

N=10

if 1:
    N=10
    x,y=np.mgrid[0:1:1/N, 0:1:1/N]
    Vx=2
    Vy=0
    A = np.sin( (2*np.pi*Vx*x)+(2*np.pi*Vy*y)+0.2*np.pi/2)
    #y = np.sin( (2*np.pi*kx*x))+np.sin((2*np.pi*ky*y))
    Ahat = np.fft.fftn(A)
    fig,axes=plt.subplots(2,2)
    ax0=axes[0][0];ax1=axes[0][1]
    ax2=axes[1][0];ax3=axes[1][1]
    ax0.imshow(A.real)
    ax1.imshow(A.imag)
    norm = mpl.colors.Normalize(vmin=0,vmax=np.abs(A).max())
    ax2.imshow(Ahat.real,norm=norm)
    ax3.imshow(Ahat.imag,norm=norm)
    fig.tight_layout()
    fig.savefig('plots_to_sort/t3')
    print('yes')
    print( np.where( np.abs(Ahat)>1e-5))

if 1:
    x,y=np.mgrid[0:1:1/N, 0:1:1/N]
    kI = np.fft.fftfreq(N)*N
    kx,ky=np.meshgrid(kI,kI)
    #yes, these are reversed from what I expect.
    ky,kx = np.meshgrid(kI,kI)

    if 0:
        Vx2=Vx+1
        Vy2=Vy+1
        Ahat = (np.abs(kx-Vx)<1e-5)*(np.abs(ky-Vy)<1e-5)*np.exp(0.4*1j*np.pi/2)
        Ahat += (np.abs(kx-Vx2)<1e-5)*(np.abs(ky-Vy2)<1e-5)*np.exp(0.3*1j*np.pi/2)
    if 1:
        r2 = kx**2+ky**2
        Ahat = np.zeros_like(r2*1j)
        Ahat[(np.abs(kx-Vx)<1e-5)*(np.abs(ky-Vy)<1e-5)]=1
    if 0:
        r2 = kx**2+ky**2
        Ahat = np.zeros_like(r2*1j)
        ok = r2>0
        Ahat[ok]=r2[ok]**-1.5
        #Ahat[ky<=0*(np.abs(kx)==0)]=0
        #Ahat[ky==0*(np.abs(kx)==0)]=0
        Ahat[ky<0]=0
        #Ahat[(np.abs(kx-0)<1e-5)*(np.abs(ky+1)<1e-5)]=1

        #Ahat[(np.abs(kx-Vx)<1e-5)*(np.abs(ky-Vy)<1e-5)]=1
    if 0:
        r2 = kx**2+ky**2
        Ahat = kx
        #Ahat[ky<0]=0
    if 0:
        phi = np.random.random(Ahat.size)
        phi.shape = Ahat.shape
        Ahatmag = np.abs(Ahat)
        Ahat = Ahatmag*np.cos(phi)+Ahatmag*np.sin(phi)*1j



    #THIS WORKS
    H = Ahat.shape[0]//2
    #Ahat[H:,H:]=np.conj(Ahat[H:0:-1,H:0:-1])
    #Ahat[:H,H:]=np.conj(Ahat[:H-1:-1,H-1::-1])
    Ahat[1:,H:] = np.conj(Ahat[:0:-1,H:0:-1])
    Ahat[0,H:] = np.conj(Ahat[0,H:0:-1])
    #Ahat[H:,0] = np.conj(Ahat[H:0:-1,0])
    #############
    print('maybe')
    print( np.where( np.abs(Ahat) > 1e-5))
    A = np.fft.ifftn(Ahat)


    fig,axes=plt.subplots(2,2)
    ax0=axes[0][0];ax1=axes[0][1]
    ax2=axes[1][0];ax3=axes[1][1]
    ax0.imshow(A.real)
    ax1.imshow(A.imag)
    print(np.abs(A.imag).sum()/np.abs(A).sum())
    #Ahat = np.roll(Ahat, Ahat.shape[0]//2,axis=0)
    #Ahat = np.roll(Ahat, Ahat.shape[0]//2,axis=1)
    ax2.imshow(Ahat.real)
    ax3.imshow(Ahat.imag)
    fig.tight_layout()
    fig.savefig('plots_to_sort/t5')

