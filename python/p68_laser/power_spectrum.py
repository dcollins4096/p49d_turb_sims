from dtools.starter1 import *
class powerspectrum():
    def __init__(self,arr):
        Nhat=np.fft.fftn(arr)
        rhohat = np.abs(Nhat)**2
        kx = np.fft.fftfreq(rhohat.shape[0])
        kabs = np.sort(np.unique(np.abs(kx)))
        kkx,kky=np.meshgrid(kx,kx)
        k = np.sqrt(kkx**2+kky**2)
        power, bins, counts =scipy.stats.binned_statistic(k.flatten(), rhohat.flatten(), bins=kabs,statistic='sum')
        bc = 0.5*(bins[1:]+bins[:-1])
        self.Nhat=Nhat  
        self.rho=arr
        self.rhohat=rhohat
        self.k = k
        self.power=power
        self.kcen=bc
