from dtools.starter1 import *
import regions
reload(regions)
import physical_values as phys
reload(phys)

import horizontal_distance as horz
reload(horz)
import equal_probability_binner as ep

class shot():
    def __init__(self,name, lines=[200,400]):
        self.name=name
        self.rho = phys.image_to_density_take2(name)
        sl = slice(*lines)
        self.std = self.rho[sl,:].std(axis=0)
        self.rhobar = self.rho[sl,:].mean(axis=0)
        self.x = phys.get_x(name)

class device():
    def __init__(self,base, lines=[200,400]):
        self.name1='%s_t1'%base
        self.name2='%s_t2'%base
        self.shot1 = shot(self.name1,lines)
        self.shot2 = shot(self.name2,lines)

    def image_device(fname):
        fig,ax=plt.subplots()

    def compute_velocity(self,rng, fname=None, nbins=16):
        print('do shock velocity')

        ok = slice(rng[0],rng[1])
        xa = self.x_1[ok]
        ya = self.rhobar_1[ok]
        xb = self.x_2[ok]
        yb = self.rhobar_2[ok]
        #horz.try2(ya,yb,method=1,fname='t1')
        if fname is not None:
            horz.try2(yb,ya,method=2,fname=fname)
        #I1 = nar(horz.ho(yb,ya))
        #dx1 = I1[:,1]-I1[:,0]
        #vel1 = phys.pixel_to_velocity(dx1)
        I2 = horz.ho2(ya=yb,yb=ya)
        dx2 = I2[:,1]-I2[:,0]
        self.vel_dist = phys.pixel_to_velocity(dx2)
        hist, cen = ep.equal_prob( self.vel_dist, nbins)
        self.vel = cen[ np.argmax(hist)]

    def bumper(self, rng,fix_shift_x=None,do_plot=False, fname='bumper.pdf'):
        from scipy.interpolate import CubicSpline
        from scipy.optimize import curve_fit
        sl = slice(*rng)
        #take the units off for the fitter.
        x_units = self.shot1.x[sl].units
        rho_units = self.shot1.rhobar.units
        x_hold = self.shot1.x[sl].v
        y_hold = self.shot1.rhobar[sl].v
        x_move = self.shot2.x[sl].v
        y_move = self.shot2.rhobar[sl].v
        interpolator = CubicSpline(x_move,y_move)
        def test_func(x,dx,dy):
            if fix_shift_x is not None:
                my_dx = fix_shift_x
            else:
                my_dx = dx
            return interpolator(x+my_dx) + dy
        self.test_func=test_func
        popt,pcov = curve_fit(test_func,x_hold,y_hold,p0=[0,0])

        self.shift_x, self.shift_rho = popt
        self.rhobar_2 = np.interp( self.shot2.x.v+popt[0], self.shot2.x.v, self.shot2.rhobar) + popt[1]
        self.rhobar_2 *= rho_units
        self.x_2 = self.shot2.x + self.shift_x*x_units
        self.rhobar_1 = self.shot1.rhobar
        self.x_1 = self.shot1.x

        if do_plot:
            fig,axes=plt.subplots(1,2)
            ax0=axes[0];ax1=axes[1]
            ax0.plot(x_hold,y_hold)
            ax0.plot(x_move,y_move)
            one_zone = x_hold[1]-x_hold[2]
            print('one zone', one_zone, 'shift',popt[0],'shift in zones', popt[0]/one_zone)
            ax0.plot( x_hold, test_func(x_hold,popt[0],popt[1]))

            ax1.plot(self.shot1.x,self.shot1.rhobar)
            ax1.plot(self.shot2.x,self.shot2.rhobar)
            ax1.plot(self.shot2.x.v, self.rhobar_2)



            fig.tight_layout()
            fig.savefig('plots_to_sort/%s'%fname)

    def atwood(self, mean_density, fname=None, gamma=5./3):

        sl = slice(mean_density[0],mean_density[1])
        rhosl = self.rhobar_1[sl]
        index = np.argmin(rhosl)
        mean_rho = self.rhobar_1[sl][index]
        sl2 = slice( mean_density[0]+index, mean_density[2])
        peak_rho = np.max( self.rhobar_1[sl2])

        print(self.vel)
        R = peak_rho/mean_rho
        self.cs = np.sqrt( gamma*self.vel**2*(1/R)*(1-1/R))
        print("cs = ",self.cs)

        if fname is not None:
            fig,axes=plt.subplots(1,3)
            ax0=axes[0]; ax1=axes[1]; ax2=axes[2]
            ax0.plot(self.rhobar_1)
            ax0.axvline(mean_density[0])
            ax0.axvline(mean_density[1])
            ax0.axvline(mean_density[2])
            ax1.plot( self.x_1[sl], self.rhobar_1[sl])
            ax1.axvline(self.x_1[sl][index])
            ax1.axhline(mean_rho)
            ax2.plot( self.rhobar_1[sl2])
            fig.tight_layout()
            fig.savefig('plots_to_sort/%s'%fname)
