from GL import *
import simulation
from scipy.optimize import curve_fit

plotdir =  dl.plotdir

def bifit_3(x, a, b, c):
    #x0 is Ms x1 is Ma
    return a+ b*x[0] + c*x[1] 
def bifit_4(x, a, b, c,d):
    #x0 is Ms x1 is Ma
    return a+ b*x[0] + c*x[1] +d*x[0]*x[1]


class beefitter():
    def __init__(self,fieldname,ms,ma,myval,fit34=3):
        self.fieldname=fieldname
        self.Ms = ms.flatten()
        self.Ma = ma.flatten()
        self.MsMa = np.row_stack([ms.flatten(),ma.flatten()])
        self.Q = myval.flatten()
        if fit34==3:
            self.fitter = bifit_3
            p0=[0,1,1]
        else:
            self.fitter = bifit_4
            p0=[0,1,1,1]
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

def plot_herd(herdlist, truth=None,predict=None):
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    verts=dt.extents()
    for nh,this_herd in enumerate(herdlist):
        ax.scatter( this_herd.Ms, this_herd.Ma, this_herd.Q)

        x1=np.linspace(this_herd.Ms.min(),this_herd.Ms.max(),10)
        y1=np.linspace(this_herd.Ma.min(),this_herd.Ma.max(),10)
        X2,Y2 = np.meshgrid(x1,y1)

        a,b,c=this_herd.Params

        surface = this_herd.Params[0]+this_herd.Params[1]*X2+this_herd.Params[2]*Y2
        ax.plot_surface(X2,Y2,surface, alpha=0.5)
        z = a + b*x1+c*y1
        verts(z)
        ax.plot(x1,y1,z)

        if truth is not None and False:
            this_truth=truth[nh]
            zt = np.zeros_like(X2)+this_truth
            ax.plot_surface(X2,Y2,zt,alpha=0.5)
            yline = (this_truth - a-b*x1)/c
            ax.plot(x1,yline,zs=this_truth)
    if predict is not None:
        z = np.linspace(verts.minmax[0],verts.minmax[1],10)
        nz=len(z)
        x = [predict[0]]*nz
        y = [predict[1]]*nz
        ax.plot(x,y,z,c='k')

    ax.set(xlabel='Ms',ylabel='Ma',zlabel='Q')
    if 1:
        # Rotate the axes and update
        angles=range(0, 360 + 1,10)
        #angles=[250]
        #angles=[0]
        for angle in angles:
            print(angle)
            # Normalize the angle to the range [-180, 180] for display
            angle_norm = (angle + 180) % 360 - 180

            # Cycle through a full rotation of elevation, then azimuth, roll, and all
            if 0:
                elev = azim = roll = 0
                if angle <= 360:
                    elev = angle_norm
                elif angle <= 360*2:
                    azim = angle_norm
                elif angle <= 360*3:
                    roll = angle_norm
                else:
                    elev = azim = roll = angle_norm
            elev = 10
            azim = angle_norm

            # Update the axis view and title
            ax.view_init(elev, azim)
            ax.set_title("elev %d axim %d"%(elev,azim))
            fig.savefig('%s/bilinear_%04d'%(dl.plotdir,angle))


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
        if len(bf1.Params)==3:
            b,c,d = bf1.Params
            g,h,k = bf2.Params
        else:
            b,c,d,WWW = bf1.Params
            g,h,k,GGG = bf2.Params
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
    fptr.write( r' Spectra & $a$ & $b$ & $c$ & $c/b$ \\'+newline)

    latex_symb = {'ClEEy_s':r'$\alpha_{EE}$', 
                  'ClBBy_s':r'$\alpha_{BB}$', 
                  'ClTTy_s':r'$\alpha_{TT}$', 
                  'density_s':r'$\alpha_\rho$', 
                  'velocity_s':r'$\alpha_v$', 
                  'magnetic_s':r'$\alpha_H$',
                  'ClEEy_a':r'$\ln A_{EE}$', 
                  'ClBBy_a':r'$\ln A_{BB}$', 
                  'ClTTy_a':r'$\ln A_{TT}$', 
                  'density_a':r'$\ln A_\rho$', 
                  'velocity_a':r'$\ln A_v$', 
                  'magnetic_a':r'$\ln A_H$'}
    for nf,field in enumerate(field_list):
        if field not in herd:
            print("Not in the herd:", field)
            continue
        line=''
        line += latex_symb.get(field,field)
        line += r'& %0.2f'%(herd[field].Params[0])
        nrm=1
        #nrm=np.abs(herd[field].Params[0]
        line += r'& %0.2f'%(herd[field].Params[1]/nrm)
        line += r'&  %0.2f'%(herd[field].Params[2]/nrm)
        line += r'&  %0.2f'%(herd[field].Params[2]/herd[field].Params[1])
        line += r'\\'
        fptr.write(line+newline)


    fptr.write( r'\hline'+newline)
    fptr.write( r'\hline'+newline)
    fptr.write( r'\end{tabular}'+newline)
    fptr.write( r'\input{table_caption}'+newline)
    fptr.write( r'\label{tab:multifit} \end{table}')
    fptr.close()

