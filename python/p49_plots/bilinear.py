from GL import *
import queb3
reload(queb3)
import davetools as dt
import sim_colors
reload(sim_colors)
import get_all_quantities as gaq
reload(gaq)
verbose=False
from collections import defaultdict
all_slopes=defaultdict(list)
from scipy.optimize import curve_fit

#read_stuff reads spectra and average quantities
#Reload to re-read, but don't do that every time.
import read_stuff as rs
spectra_dict = rs.spectra_dict
import read_avg_quan as raq
quan3=raq.quan3

plotdir =  "/home/dccollins/PigPen"

#def bifit(x, a, b, c, d):
def bifit_4(x, a, b, c,d):
    #x0 is Ms x1 is Ma
    return a+ b*x[0] + c*x[1] + d*x[0]*x[1]
    #return a+ b*x[0] + c*x[1] 
def bifit_3(x, a, b, c):
    #x0 is Ms x1 is Ma
    return a+ b*x[0] + c*x[1] 
def bifit_2(x,  b, c):
    #x0 is Ms x1 is Ma
    return b*x[0] + c*x[1] 
def do_bifit(xdata, ydata, p0):
    fitParams, fitCovariances = curve_fit(bifit, xdata, ydata, p0)
    #print(' fit coefficients:\n', fitParams)
    return fitParams, fitCovariances
class beefitter():
    def __init__(self,fieldname,ms,ma,myval):
        self.fieldname=fieldname
        self.MsMa = np.row_stack([ms.flatten(),ma.flatten()])
        self.Q = myval.flatten()
        p0=[0,1,1]
        self.fitter = bifit_3
        #p0=[1,1]
        #self.fitter = bifit_2
        self.Params, self.Cov = curve_fit( self.fitter, self.MsMa, self.Q, p0)
    def plot(self):
        plt.clf()
        Qprime = self.fitter(self.MsMa, *self.Params)
        for n in range( len( self.Q)):
            plt.scatter( self.MsMa[0][n], self.Q[n], marker=list(sim_colors.markerlist)[n], c=sim_colors.colorlist[n])
            plt.scatter( self.MsMa[0][n], Qprime[n], marker=list(sim_colors.markerlist)[n], c=sim_colors.colorlist[n])

        plt.savefig('%s/linear_%s.png'%(plotdir,self.fieldname))

def do_morebifit(fieldname, ms, ma, myval):
    x=np.row_stack([ms.flatten(),ma.flatten()])
    y=myval.flatten()
    p0 = [0,1,1,0]
    print(fieldname)
    fitParams, fitCovariances= do_bifit(x,y, p0)
    f = fitParams
    cov = fitCovariances
#    ys = f[0] + f[1]*myx + f[2]*myy+f[3]*myx*myy
    Fptr = h5py.File("bi_fit_params.h5",'a')
    Fptr.create_dataset(fieldname, data=f)
    Fptr.close()
    Fptr = h5py.File("bi_fit_covariances.h5",'a')
    Fptr.create_dataset(fieldname+"_cov", data=cov)
    Fptr.close()


def pairwise( bf1, bf2, truth1, truth2):
    if 0:
        p = bf1.Params
        q = bf2.Params
        Ms =  (q[2]*(p[0]-truth1)+p[2]*(q[0]-truth2))/(p[2]*q[1]-p[1]*q[2])
        Ma = -(q[1]*(p[0]-truth1)+p[1]*(q[0]-truth2))/(p[2]*q[1]-p[1]*q[2])
    else:
        # a = b + c x + d y 
        # f = g + h x + k y
        # solve pair wise
        a = truth1
        f = truth2
        b,c,d = bf1.Params
        g,h,k = bf2.Params
        x = (-a*k+b*k+d*f-d*g)/(d*h-c*k)
        y = ( h*(b-a)+c*(f-g))/(c*k-d*h)
        Ms=x;Ma=y
    return Ms, Ma

