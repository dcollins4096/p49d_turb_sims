from dtools.starter1 import *
import regions
reload(regions)
import physical_values as phys
reload(phys)

import horizontal_distance as horz
reload(horz)


class shot():
    def __init__(self,name, lines=[200,400]):
        self.name=name
        self.rho = phys.useful_values_take2(name)
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

    def compute_velocity(self,rng):
        print('do shock velocity')

        ok = slice(rng[0],rng[1])
        xa = self.x_1[ok]
        ya = self.rhobar_1[ok]
        xb = self.x_2[ok]
        yb = self.rhobar_2[ok]
        #horz.try2(ya,yb,method=1,fname='t1')
        #horz.try2(yb,ya,method=2,fname='t2')
        #I1 = nar(horz.ho(yb,ya))
        #dx1 = I1[:,1]-I1[:,0]
        #vel1 = phys.pixel_to_velocity(dx1)
        I2 = horz.ho2(ya=yb,yb=ya)
        dx2 = I2[:,1]-I2[:,0]
        self.vel = phys.pixel_to_velocity(dx2)





    def bumper(self, rng,do_plot=False):
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
            return interpolator(x+dx) + dy
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
            ax1.plot(self.shot2.x.v, interp)



            fig.tight_layout()
            fig.savefig('plots_to_sort/bumper.pdf')

