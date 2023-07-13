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
        self.Ms = ms.flatten()
        self.Ma = ma.flatten()
        self.MsMa = np.row_stack([ms.flatten(),ma.flatten()])
        self.Q = myval.flatten()
        p0=[0,1,1]
        self.fitter = bifit_3
        #p0=[1,1]
        #self.fitter = bifit_2
        self.Params, self.Cov = curve_fit( self.fitter, self.MsMa, self.Q, p0)
    def plot1(self):
        plt.clf()
        Qprime = self.fitter(self.MsMa, *self.Params)
        for n in range( len( self.Q)):
            plt.scatter( self.MsMa[0][n], self.Q[n], marker=list(sim_colors.markerlist)[n], c=sim_colors.colorlist[n])
            plt.scatter( self.MsMa[0][n], Qprime[n], marker=list(sim_colors.markerlist)[n], c=sim_colors.colorlist[n])
        plt.savefig('%s/linear_%s.png'%(plotdir,self.fieldname))
    def plot2(self):
        plt.clf()
        Qprime = self.fitter(self.MsMa, *self.Params)
        for n in range( len( self.Q)):
            plt.scatter( self.MsMa[0][n], self.Q[n], marker=list(sim_colors.markerlist)[n], c=sim_colors.colorlist[n])
            #plt.scatter( self.MsMa[0][n], Qprime[n], marker=list(sim_colors.markerlist)[n], c=sim_colors.colorlist[n])
        for ma in self.Ma:
            tmp_ms = nar([0]+list(self.Ms))
            this_ms_ma = [tmp_ms,ma]
            #plt.plot( tmp_ms, self.fitter(this_ms_ma, *self.Params))
            print(self.Params)


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

def write_tex(herd, fname, field_list=None):
    if field_list is None:
        field_list = herd.keys()
    fptr=open(fname,'w')
    newline='\n'
    fptr.write( r'\begin{table}'+newline)
    fptr.write( r'\begin{tabular}{lrrrr}'+newline)
    fptr.write( r'\hline'+newline)
    #fptr.write( r' Spectra & $a$ & $b/\left|a\right|$ & $c/\left|a\right|$ & $d/\left|a\right|$ \\'+newline)
    #fptr.write( r' Spectra & $a$ & $b/\left|a\right|$ & $c/\left|a\right|$ \\'+newline)
    fptr.write( r' Spectra & $a$ & $b$ & $c$ \\'+newline)

    latex_symb = {'avg_clees':r'$\alpha_{EE}$', 'avg_clbbs':r'$\alpha_{BB}$', 
                  'avg_cltts':r'$\alpha_{TT}$', 'avg_ds':r'$\alpha_\rho$', 
                  'avg_vs':r'$\alpha_v$', 'avg_hs':r'$\alpha_H$',
                  'avg_cleea':r'$\ln A_{EE}$', 'avg_clbba':r'$\ln A_{BB}$', 
                  'avg_cltta':r'$\ln A_{TT}$', 'avg_da':r'$\ln A_\rho$', 
                  'avg_va':r'$\ln A_v$', 'avg_ha':r'$\ln A_H$'}
    for nf,field in enumerate(field_list):
        if field not in herd:
            print("Not in the herd:", field)
            continue
        line=''
        line += latex_symb[field]
        line += r'& %0.2f'%(herd[field].Params[0])
        nrm=1
        #nrm=np.abs(herd[field].Params[0]
        line += r'& %0.2f'%(herd[field].Params[1]/nrm)
        line += r'&  %0.2f'%(herd[field].Params[2]/nrm)
        line += r'\\'
        fptr.write(line+newline)


    fptr.write( r'\hline'+newline)
    fptr.write( r'\hline'+newline)
    fptr.write( r'\end{tabular}'+newline)
    fptr.write( r'\input{table_caption}'+newline)
    fptr.write( r'\label{tab:multifit} \end{table}')
    fptr.close()

