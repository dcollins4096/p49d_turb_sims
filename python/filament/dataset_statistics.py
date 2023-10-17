import yt
import os
import time
import hessian
import histogram
import numpy as np
import matplotlib.pyplot as plt
import pcolormesh_helper as pch

def angle_weigh(angle):
       res = np.sin(angle)
       res[res == 0] = np.inf
       res = 1 / res
       return res

def write_number(num, filename):
        f = open(filename, "w")
        f.write("%.15f"%num)
        f.close()

def write_histogram(hist, filename):
        f = open(filename, "w")
        print("Writing to file %s..."%filename)
        for entry in hist:
                line = ' '.join((str(e) for e in entry))
                f.write(line)
                f.write("\n")
        f.close()
        print("done!\n")

def write_1D(arr, filename):
        f = open(filename, "w")
        print("Writing to file %s..."%filename)
        for entry in arr:
                f.write(entry)
                f.write("\n")
        f.close()
        print("done!\n")

def write_2D(arr1, arr2, filename):
        f = open(filename, "w")
        print("Writing to file %s..."%filename)
        for i in range(arr1.size):
                f.write("%.15f %.15f\n"%(arr1[i],arr2[i]))
        f.close()
        print("done!\n")

def write_fourier2D(freq_x, freq_y, fourier, filename):
        f = open(filename, "w")
        print("Writing to file %s..."%filename)
        for i in range(fourier.shape[0]):
               for j in range(fourier.shape[1]):
                      f.write("%.15f %.15f %.15f\n"%(freq_x[i], freq_y[j], fourier[i,j]))
        f.close()
        print("done!\n")

def blob_statistics(eigenvalues, density):
       mask = eigenvalues[...,2] < 0
       weight = eigenvalues[...,2] / eigenvalues[...,0]
       total = np.prod(eigenvalues.shape)
       count = np.sum(mask, axis = None)
       #print(count / total)
       count = np.sum(mask * weight, axis = None)
       #print(count / total)
       plt.clf()
       plt.imshow((density * mask).sum(axis=0))
       plt.savefig("blobs.png")
       
       plt.clf()
       plt.imshow((density * mask * weight).sum(axis=0))
       plt.savefig("blobs_weighted.png")

def cos2acos(x):
       return 2*x**2 - 1

def sin2acos(x):
       return 2*x*np.sqrt(1 - x**2)

def getQU(density, Bx, By, Bz, ax = "x"):
       # gamma is the angle between B and the POS. gamma2 is the angle between B and the LOS
       # cos(gamma)^2 = sin(pi/2 - gamma)^2 = sin^2(gamma2) = 1 - cos^2(gamma2) = 1 - (BLOS / |B|)^2
       # psi: cos(psi) = By / |BPOS|, sin(psi) = Bz / |BPOS|
       # cos(2*arccos(x)) = 2 x^2 - 1
       # sin(2*arccos(x)) = 2 x sqrt(1 - x^2)
       magB2 = Bx**2 + By**2 + Bz**2
       if ax == "x":
              BPOSmag = np.sqrt(By**2 + Bz**2)
              # Bz dominated: q = 1, u = 0
              # By dominated: q = -1, u = 0
              # By and Bz are equal: q = 0, u = 1
              #q = (density * cos2acos(Bz / BPOSmag) * (1 - (Bx / magB)**2)).sum(axis=0) # up / right
              #u = (density * sin2acos(By / BPOSmag) * (1 - (Bx / magB)**2)).sum(axis=0) # diagonal
              q = (density * (By**2 - Bz**2) / magB2).sum(axis=0)
              u = 2 * (density * (By * Bz) / magB2).sum(axis=0)
       elif ax == "y":
              BPOSmag = np.sqrt(Bx**2 + Bz**2)
              #q = (density * cos2acos(Bx / BPOSmag) * (1 - (By / magB)**2)).sum(axis=1)
              #u = (density * sin2acos(Bz / BPOSmag) * (1 - (By / magB)**2)).sum(axis=1)
              q = (density * (Bx**2 - Bz**2) / magB2).sum(axis=1)
              u = 2 * (density * (Bx * Bz) / magB2).sum(axis=1)
       else:
              BPOSmag = np.sqrt(Bx**2 + By**2)
              # rewrite these
              #q = (density * cos2acos(Bx / BPOSmag) * (1 - (Bz / magB)**2)).sum(axis=2)
              #u = (density * sin2acos(By / BPOSmag) * (1 - (Bz / magB)**2)).sum(axis=2) # 2 Bx |By| / magB^2
              q = (density * (Bx**2 - By**2) / magB2).sum(axis=2)
              u = 2 * (density * (Bx * By) / magB2).sum(axis=2)
       # 0.5 * atan2(u = 0, q = +1) = 0
       # 0.5 * atan2(u = 0, q = -1) = 90
       # 0.5 * atan2(u = +1, q = 0) = 45
       # 0.5 * atan2(u = -1, q = 0) = -45
       psi = 0.5 * np.arctan2(u, q) #atan(y/x)
       plt.clf()
       plt.imshow(q)
       plt.savefig(pics + "Q_" + ax + ".png")
       plt.clf()
       plt.imshow(u)
       plt.savefig(pics + "U_" + ax + ".png")
       plt.clf()
       plt.imshow(psi)
       plt.savefig(pics + "psi_" + ax + ".png")
       return psi, q, u

def getEB(q, u, ax, postfix):
       n = q.size
       nx = q.shape[0]
       ny = q.shape[1]
       fourierQ = np.fft.fft2(q) / n
       fourierU = np.fft.fft2(u) / n
       kx = 2 * np.pi * np.fft.fftfreq(q.shape[0])
       ky = 2 * np.pi * np.fft.fftfreq(q.shape[1])
       kx_arr, ky_arr = np.meshgrid(kx ,ky)
       #fourierAngle = np.arctan2(ky_arr, kx_arr) # atan(ky/kx)
       # (kx + i*ky) / |k|
       # cos(angle) = kx / |k|
       # sin(angle) = ky / |k|
       # cos(2angle) = cos(angle)^2 - sin(angle)^2 = (kx^2 - ky^2) / |k|^2
       # sin(2angle) = 2 sin(angle) cos(angle) = 2 kx ky / |k|^2
       kmag2 = kx_arr**2 + ky_arr**2
       c = (kx_arr**2 - ky_arr**2) / kmag2
       s = 2 * kx_arr * ky_arr / kmag2
       modeEfourier = fourierQ * c + fourierU * s
       modeBfourier = fourierU * c - fourierQ * s
       modeEfourier[0,0] = fourierQ[0,0]
       modeBfourier[0,0] = fourierU[0,0]
       modeE = np.fft.ifft2(modeEfourier)
       modeB = np.fft.ifft2(modeBfourier)
       plt.clf()
       plt.imshow(np.real(modeE))
       plt.savefig(pics + "E_" + ax + postfix + ".png")
       plt.clf()
       plt.imshow(np.real(modeB))
       plt.savefig(pics + "B_" + ax + postfix + ".png")
       #write_fourier2D(kx, ky, np.real(modeEfourier), stats + "mode_re_E_fourier_" + ax + ".txt")
       #write_fourier2D(kx, ky, np.imag(modeEfourier), stats + "mode_im_E_fourier_" + ax + ".txt")
       #write_fourier2D(kx, ky, np.real(modeBfourier), stats + "mode_re_B_fourier_" + ax + ".txt")
       #write_fourier2D(kx, ky, np.imag(modeBfourier), stats + "mode_im_B_fourier_" + ax + ".txt")
       return kx, ky, modeEfourier, modeBfourier, modeE, modeB

def fermi(energy, temp):
       return 1 / (1 + np.exp(energy / temp))

def fermi2(energy, temp):
       return 1 / (1 + np.exp((energy - temp**2) / temp))

def smoother(energy, temp):
       poz = energy > 0
       return np.exp( - (energy / temp)**2) * poz + 1 - poz

def mask3D(eigenvalues):
       return (eigenvalues[...,1] < 0) * (abs(eigenvalues[...,0]) > abs(eigenvalues[...,2])) * (abs(eigenvalues[...,1]) > abs(eigenvalues[...,2]))

def mask3D_smooth(eigenvalues, temp):
       return fermi(eigenvalues[...,1], temp) * fermi(abs(eigenvalues[...,2]) - abs(eigenvalues[...,0]), temp) * fermi(abs(eigenvalues[...,2]) - abs(eigenvalues[...,1]), temp)

def mask3D_smooth2(eigenvalues, temp):
       return fermi2(eigenvalues[...,1], temp) * fermi2(abs(eigenvalues[...,2]) - abs(eigenvalues[...,0]), temp) * fermi2(abs(eigenvalues[...,2]) - abs(eigenvalues[...,1]), temp)

def mask3D_smooth3(eigenvalues, temp):
       return smoother(eigenvalues[...,1], temp) * smoother(abs(eigenvalues[...,2]) - abs(eigenvalues[...,0]), temp) * smoother(abs(eigenvalues[...,2]) - abs(eigenvalues[...,1]), temp)

def filament_statistics_3D(eigenvalues, eigenvectors, density, Bx, By, Bz):
       mask = mask3D(eigenvalues)
       weight = eigenvalues[...,1] / eigenvalues[...,0]
       weight2 = np.sqrt(abs(eigenvalues[...,0] * eigenvalues[...,1]))[mask]
       total = np.prod(eigenvalues.shape)
       count = np.sum(mask, axis = None)
       #(count / total)
       count = np.sum(mask * weight, axis = None)
       #print(count / total)
       magB = np.sqrt(Bx**2 + By**2 + Bz**2)

       dot_prod = Bx * eigenvectors[...,2,0] + By * eigenvectors[...,2,1] + Bz * eigenvectors[...,2,2]
       angle = np.arccos(abs(dot_prod) / magB)[mask]

       #pch.simple_phase(np.log(density).ravel(), angle.ravel())
       #pch.simple_phase(np.log(density[mask]), angle[mask])

       density_masked = density[mask]
       angleW = angle_weigh(angle)

       hist = histogram.histogram2D_weighted(np.log(density_masked), angle, angleW, 40, 90)
       write_histogram(hist, stats + "logrho_angle_3D_masked.txt")

       hist = histogram.histogram2D_weighted(np.log(density_masked), angle, angleW * density_masked, 40, 90)
       write_histogram(hist, stats + "logrho_angle_3D_masked_weighted_M.txt")

       hist = histogram.histogram2D_weighted(np.log(density_masked), angle, angleW * weight2, 40, 90)
       write_histogram(hist, stats + "logrho_angle_3D_masked_weighted_lambda.txt")

       hist = histogram.histogram2D_weighted(np.log(density_masked), angle, angleW * density_masked * weight2, 40, 90)
       write_histogram(hist, stats + "logrho_angle_3D_masked_weighted_lambda_M.txt")
       '''
       hist = histogram.histogram_weighted(angle, angleW, 0.01)
       write_histogram(hist, stats + "filaments_angle_3D.txt")
       hist = histogram.histogram_weighted(angle, angleW * density_masked, 0.01)
       write_histogram(hist, stats + "filaments_angle_3D_weighted_M.txt")
       hist = histogram.histogram_weighted(angle, weight2 * angleW, 0.01)
       write_histogram(hist, stats + "filaments_angle_3D_weighted_lambda.txt")
       hist = histogram.histogram_weighted(angle, weight2 * angleW * density_masked, 0.01)
       write_histogram(hist, stats + "filaments_angle_3D_weighted_lambda_M.txt")

       sigma = np.sum(angleW * angle**2, axis = None)
       denom = np.sum(angleW, axis = None)
       write_number(np.sqrt(sigma / denom), stats + "filament_angle_3D_sigma.txt")

       plt.clf()
       plt.imshow(density.sum(axis=0))
       plt.savefig(pics + "density_projection_x.png")

       plt.clf()
       plt.imshow(density.sum(axis=1))
       plt.savefig(pics + "density_projection_y.png")

       plt.clf()
       plt.imshow(density.sum(axis=2))
       plt.savefig(pics + "density_projection_z.png")

       #print("Filaments")
       plt.clf()
       plt.imshow((density * mask).sum(axis=0))
       plt.savefig(pics + "filaments_3D_x.png")

       plt.clf()
       plt.imshow((density * mask).sum(axis=1))
       plt.savefig(pics + "filaments_3D_y.png")

       plt.clf()
       plt.imshow((density * mask).sum(axis=2))
       plt.savefig(pics + "filaments_3D_z.png")

       plt.clf()
       plt.imshow((density * mask * weight).sum(axis=0))
       plt.savefig(pics + "filaments_3D_x_weighted.png")

       plt.clf()
       plt.imshow((density * mask * weight).sum(axis=1))
       plt.savefig(pics + "filaments_3D_y_weighted.png")

       plt.clf()
       plt.imshow((density * mask * weight).sum(axis=2))
       plt.savefig(pics + "filaments_3D_z_weighted.png")

       weight2 = np.sqrt(abs(eigenvalues[...,0] * eigenvalues[...,1]))
       plt.clf()
       plt.imshow((density * mask * weight2).sum(axis=0))
       plt.savefig(pics + "filaments_3D_x_weighted2.png")

       plt.clf()
       plt.imshow((density * mask * weight2).sum(axis=1))
       plt.savefig(pics + "filaments_3D_y_weighted2.png")

       plt.clf()
       plt.imshow((density * mask * weight2).sum(axis=2))
       plt.savefig(pics + "filaments_3D_z_weighted2.png")
       '''

def angle_diff(angle1, angle2): # unsigned angle difference
       angle1[angle1>np.pi/2] = -np.pi + angle1[angle1>np.pi/2]
       angle1[angle1<-np.pi/2] = np.pi + angle1[angle1<-np.pi/2]
       
       angle2[angle2>np.pi/2] = -np.pi + angle2[angle2>np.pi/2]
       angle2[angle2<-np.pi/2] = np.pi + angle2[angle2<-np.pi/2]

       diff = abs(angle1 - angle2)
       diff[diff>np.pi/2] = np.pi - diff[diff>np.pi/2]

       return diff

def filament_statistics_2D(eigenvalues, eigenvectors, density, Bx, By, Bz, axis):
       mask = (eigenvalues[...,0] < 0) * (abs(eigenvalues[...,0]) > abs(eigenvalues[...,1]))
       weight = eigenvalues[...,1] / eigenvalues[...,0]
       weight2 = abs(eigenvalues[...,0])
       total = np.prod(eigenvalues.shape)
       count = np.sum(mask, axis = None)
       #print(count / total)
       count = np.sum(mask * weight2, axis = None)
       #print(count / total)
       
       if axis == "x":
              density_projection = density.sum(axis=0)
       elif axis == "y":
              density_projection = density.sum(axis=1)
       else:
              density_projection = density.sum(axis=2)
       
       #psi_B, q, u = getQU(density, Bx, By, Bz, ax=axis)
       #getEB(q, u, axis)

       # 1D histograms of Q/U angle
       hist = histogram.histogram(psi_B[mask], 0.025)
       write_histogram(hist, stats + "filaments_angle_QU_2D_" + axis + ".txt")
       hist = histogram.histogram_weighted(psi_B[mask], density_projection[mask], 0.025)
       write_histogram(hist, stats + "filaments_angle_QU_2D_" + axis + "_weighted_M.txt")
       hist = histogram.histogram_weighted(psi_B[mask], weight2[mask], 0.025)
       write_histogram(hist, stats + "filaments_angle_QU_2D_" + axis + "_weighted_lambda.txt")
       hist = histogram.histogram_weighted(psi_B[mask], weight2[mask] * density_projection[mask], 0.025)
       write_histogram(hist, stats + "filaments_angle_QU_2D_" + axis + "_weighted_lambda_M.txt")
       plt.clf()
       plt.imshow(psi_B)
       plt.savefig(pics + "psi_B" + axis + ".png")
       
       # 2D histograms of logrho and Q/U angle
       hist = histogram.histogram2D(np.log(density_projection[mask]), psi_B[mask], 40, 60)
       write_histogram(hist, stats + "logrho_angle_QU_2D_" + axis + ".txt")
       hist = histogram.histogram2D_weighted(np.log(density_projection[mask]), psi_B[mask], density_projection[mask], 40, 60)
       write_histogram(hist, stats + "logrho_angle_QU_2D_" + axis + "_weighted_M.txt")
       hist = histogram.histogram2D_weighted(np.log(density_projection[mask]), psi_B[mask], weight2[mask], 40, 60)
       write_histogram(hist, stats + "logrho_angle_QU_2D_" + axis + "_weighted_lambda.txt")
       hist = histogram.histogram2D_weighted(np.log(density_projection[mask]), psi_B[mask], weight2[mask] * density_projection[mask], 40, 60)
       write_histogram(hist, stats + "logrho_angle_QU_2D_" + axis + "_weighted_lambda_M.txt")
       
       # atan2(0,1) = 0
       # atan2(0,-1) = 180
       # atan2(1,0) = 90
       # atan2(-1,0) = -90
       psi_fil = np.arctan2(eigenvectors[...,1,1], eigenvectors[...,1,0]) # v pointing right is angle = 0
       psi_fil[psi_fil>np.pi/2] = -np.pi + psi_fil[psi_fil>np.pi/2]
       psi_fil[psi_fil<-np.pi/2] = np.pi + psi_fil[psi_fil<-np.pi/2]

       # 1D histograms of filament angle
       hist = histogram.histogram(psi_fil[mask], 0.025)
       write_histogram(hist, stats + "filaments_angle_fil_2D_" + axis + ".txt")
       hist = histogram.histogram_weighted(psi_fil[mask], density_projection[mask], 0.025)
       write_histogram(hist, stats + "filaments_angle_fil_2D_" + axis + "_weighted_M.txt")
       hist = histogram.histogram_weighted(psi_fil[mask], weight2[mask], 0.025)
       write_histogram(hist, stats + "filaments_angle_fil_2D_" + axis + "_weighted_lambda.txt")
       hist = histogram.histogram_weighted(psi_fil[mask], weight2[mask] * density_projection[mask], 0.025)
       write_histogram(hist, stats + "filaments_angle_fil_2D_" + axis + "_weighted_lambda_M.txt")
       
       # 2D histograms of logrho and filament angle
       hist = histogram.histogram2D(np.log(density_projection[mask]), psi_fil[mask], 40, 60)
       write_histogram(hist, stats + "logrho_angle_fil_2D_" + axis + ".txt")
       hist = histogram.histogram2D_weighted(np.log(density_projection[mask]), psi_fil[mask], density_projection[mask], 40, 60)
       write_histogram(hist, stats + "logrho_angle_fil_2D_" + axis + "_weighted_M.txt")
       hist = histogram.histogram2D_weighted(np.log(density_projection[mask]), psi_fil[mask], weight2[mask], 40, 60)
       write_histogram(hist, stats + "logrho_angle_fil_2D_" + axis + "_weighted_lambda.txt")
       hist = histogram.histogram2D_weighted(np.log(density_projection[mask]), psi_fil[mask], weight2[mask] * density_projection[mask], 40, 60)
       write_histogram(hist, stats + "logrho_angle_fil_2D_" + axis + "_weighted_lambda_M.txt")

       psi_diff = angle_diff(psi_B, psi_fil)
       # 1D histograms of angle difference between filaments and Q/U
       hist = histogram.histogram(psi_diff[mask], 0.025)
       write_histogram(hist, stats + "filaments_angle_2D_" + axis + ".txt")
       hist = histogram.histogram_weighted(psi_diff[mask], density_projection[mask], 0.025)
       write_histogram(hist, stats + "filaments_angle_2D_" + axis + "_weighted_M.txt")
       hist = histogram.histogram_weighted(psi_diff[mask], weight2[mask], 0.025)
       write_histogram(hist, stats + "filaments_angle_2D_" + axis + "_weighted_lambda.txt")
       hist = histogram.histogram_weighted(psi_diff[mask], weight2[mask] * density_projection[mask], 0.025)
       write_histogram(hist, stats + "filaments_angle_2D_" + axis + "_weighted_lambda_M.txt")
       
       # 2D histograms of angle difference between filaments and Q/U
       hist = histogram.histogram2D(np.log(density_projection[mask]), psi_diff[mask], 40, 60)
       write_histogram(hist, stats + "logrho_angle_2D_" + axis + ".txt")
       hist = histogram.histogram2D_weighted(np.log(density_projection[mask]), psi_diff[mask], density_projection[mask], 40, 60)
       write_histogram(hist, stats + "logrho_angle_2D_" + axis + "_weighted_M.txt")
       hist = histogram.histogram2D_weighted(np.log(density_projection[mask]), psi_diff[mask], weight2[mask], 40, 60)
       write_histogram(hist, stats + "logrho_angle_2D_" + axis + "_weighted_lambda.txt")
       hist = histogram.histogram2D_weighted(np.log(density_projection[mask]), psi_diff[mask], weight2[mask] * density_projection[mask], 40, 60)
       write_histogram(hist, stats + "logrho_angle_2D_" + axis + "_weighted_lambda_M.txt")

       #write_2D((psi_B[mask]).ravel(), (psi_fil[mask]).ravel(), stats + "angle_QU_filament_" + axis + ".txt")

       plt.clf()
       plt.imshow(density_projection * mask)
       plt.savefig(pics + "filaments_2D_" + axis + ".png")

       plt.clf()
       plt.imshow(density_projection * mask * weight)
       plt.savefig(pics + "filaments_2D_" + axis + "_weighted.png")

       plt.clf()
       plt.imshow(density_projection * mask * weight2)
       plt.savefig(pics + "filaments_2D_" + axis + "_weighted2.png")

       return density_projection * mask * weight2

def sh_weight(x): # = 1 when x = 0, = 0 when x = 1
       return 1 - x

def sheet_statistics(eigenvalues, eigenvectors, density):
       mask = (eigenvalues[...,0] < 0) * (abs(eigenvalues[...,0]) > abs(eigenvalues[...,1])) * (abs(eigenvalues[...,0]) > abs(eigenvalues[...,2]))
       weight = sh_weight(abs(eigenvalues[...,1] / eigenvalues[...,0])) * sh_weight(abs(eigenvalues[...,2] / eigenvalues[...,0]))
       total = np.prod(eigenvalues.shape)
       count = np.sum(mask, axis = None)
       #print(count / total)

       dot_prod = Bx * eigenvectors[...,0,0] + By * eigenvectors[...,0,1] + Bz * eigenvectors[...,0,2]
       magB = np.sqrt(Bx**2 + By**2 + Bz**2)
       angle = np.arccos(abs(dot_prod) / magB)

       hist = histogram.histogram(angle.ravel()[mask.ravel()], 0.0025)

       write_histogram(hist, "angle_sheets.txt")

       hist = histogram.histogram_weighted(angle.ravel()[mask.ravel()], weight.ravel()[mask.ravel()], 0.0025)

       write_histogram(hist, "angle_sheets_weighted.txt")

       plt.clf()
       plt.imshow((density * mask).sum(axis=0))
       plt.savefig("sheets.png")
       
       plt.clf()
       plt.imshow((density * mask * weight).sum(axis=0))
       plt.savefig("sheets_weighted.png")

def unit_inverse(x, y, z, M):
       rho = np.sqrt(x**2+y**2)
       r = np.sqrt(x**2+y**2+z**2)
       M[..., 0, 0] = - y/rho
       M[..., 0, 1] = x/rho
       M[..., 0, 2] = 0
       
       M[..., 1, 0] = - x * z / (rho * r)
       M[..., 1, 1] = - y * z / (rho * r)
       M[..., 1, 2] = rho / r
       
       M[..., 2, 0] = x / r
       M[..., 2, 1] = y / r
       M[..., 2, 2] = z / r

def matrix_multiply(M, x, y, z):
       # xvec = x[..., np.newaxis]
       # xvec = np.repeat(xvec,3,axis=-1)
       # xvec = xvec[...,np.newaxis]
       xvec = np.stack((x, y, z), axis=-1)[..., np.newaxis] # (..., 3, 1) ?
       #print(M.shape)
       #print(xvec.shape)
       return M @ xvec

def stats_Bv(Bx, By, Bz, vx, vy, vz, rho):
       Bx_avg = np.mean(Bx)
       By_avg = np.mean(By)
       Bz_avg = np.mean(Bz)
       magB = np.sqrt(Bx**2 + By**2 + Bz**2)
       magV = np.sqrt(vx**2 + vy**2 + vz**2)
       M = vx[..., np.newaxis, np.newaxis]
       M = np.repeat(M, 3, axis=-2)
       M = np.repeat(M, 3, axis=-1) # (512, 512, 512, 3, 3)
       unit_inverse(vx, vy, vz, M)
       Brot = matrix_multiply(M, Bx - Bx_avg, By - By_avg, Bz - Bz_avg) # (512, 512, 512, 3, 1)
       angle = np.arccos((Bx * vx + By * vy + Bz * vz) / (magB * magV))

       hist = histogram.histogram_weighted(angle.ravel(), (angle_weigh(angle)).ravel(), 0.01)
       write_histogram(hist, stats + "angle_B_v.txt")

       hist = histogram.histogram_weighted(angle.ravel(), (angle_weigh(angle) * rho).ravel(), 0.01)
       write_histogram(hist, stats + "angle_B_v_weighted_M.txt")
       
       hist = histogram.histogram2D_weighted(np.log(rho.ravel()), angle.ravel(), angle_weigh(angle.ravel()), 50, 180)
       write_histogram(hist, stats + "logrho_B_v_angle.txt")
       
       hist = histogram.histogram2D_weighted(np.log(rho.ravel()), angle.ravel(), angle_weigh(angle.ravel()) * rho.ravel(), 50, 180)
       write_histogram(hist, stats + "logrho_B_v_angle_weighted_M.txt")
       
       v2 = vx * vx + vy * vy + vz * vz

       hist = histogram.histogram2D_weighted(np.log(rho.ravel()), angle.ravel(), angle_weigh(angle.ravel()) * 0.5 * rho.ravel() * v2.ravel(), 50, 180)
       write_histogram(hist, stats + "logrho_B_v_angle_weighted_KE.txt")

       '''
       write_histogram(hist, stats + "angle_B_v.txt")

       hist = histogram.histogram(Bx.ravel() - Bx_avg, 0.0025)

       write_histogram(hist, "Bx.txt")

       hist = histogram.histogram(By.ravel() - By_avg, 0.0025)

       write_histogram(hist, "By.txt")

       hist = histogram.histogram(Bz.ravel() - Bz_avg, 0.0025)

       write_histogram(hist, "Bz.txt")

       magB = np.sqrt((Bx-Bx_avg)**2 + (By-By_avg)**2 + (Bz-Bz_avg)**2)
       hist = histogram.histogram(magB.ravel(), 0.0025)

       write_histogram(hist, "abs_B.txt")


       hist = histogram.histogram(Brot[...,0,0].ravel(), 0.0025)

       write_histogram(hist, "Bx_rel.txt")

       hist = histogram.histogram(Brot[...,1,0].ravel(), 0.0025)

       write_histogram(hist, "By_rel.txt")

       hist = histogram.histogram(Brot[...,2,0].ravel(), 0.0025)

       write_histogram(hist, "Bz_rel.txt")
       '''


def density_stats(rho, vx, vy, vz):
       max_power = 100

       #projections histograms
       n = rho.shape[0]
       rho_proj = (rho.sum(axis = 0) / n)
       logrho_proj = np.log(rho_proj)
       hist = histogram.histogram(logrho_proj.ravel(), 0.01)
       write_histogram(hist, stats + "hist1D_logrho_x.txt")
       hist = histogram.histogram_weighted(logrho_proj.ravel(), rho_proj.ravel(), 0.01)
       write_histogram(hist, stats + "hist1D_logrho_x_weighted_M.txt")
       n = rho_proj.size
       s = logrho_proj.sum() / n
       s2 = (logrho_proj**2).sum() / n
       rho2D = rho_proj.sum() / n
       write_number(s, stats + "logrho_mean_2D_x.txt")
       write_number(np.sqrt(s2 - s**2), stats + "logrho_sigma_2D_x.txt")
       write_number(rho2D, stats + "rho_2D_x.txt")
       for power in range(1,max_power + 1):
              rho_pow = (rho_proj**power).sum() / n
              write_number(rho_pow, stats + "rho_pow_" + str(power) + "_2D_x.txt")
       
       #1D projections
       n = rho.shape[0]
       rho_proj1D = (rho_proj.sum(axis = 0) / n).ravel()
       n = rho_proj1D.size
       s = np.log(rho_proj1D).sum() / n
       s2 = (np.log(rho_proj1D)**2).sum() / n
       write_number(s, stats + "logrho_mean_1D_xy.txt")
       write_number(np.sqrt(s2 - s**2), stats + "logrho_sigma_1D_xy.txt")
       for power in range(1,max_power + 1):
              rho_pow = (rho_proj1D**power).sum() / n
              write_number(rho_pow, stats + "rho_pow_" + str(power) + "_1D_xy.txt")
       
       n = rho.shape[0]
       rho_proj1D = (rho_proj.sum(axis = 1) / n).ravel()
       n = rho_proj1D.size
       s = np.log(rho_proj1D).sum() / n
       s2 = (np.log(rho_proj1D)**2).sum() / n
       write_number(s, stats + "logrho_mean_1D_xz.txt")
       write_number(np.sqrt(s2 - s**2), stats + "logrho_sigma_1D_xz.txt")
       for power in range(1,max_power + 1):
              rho_pow = (rho_proj1D**power).sum() / n
              write_number(rho_pow, stats + "rho_pow_" + str(power) + "_1D_xz.txt")
       
       n = rho.shape[1]
       rho_proj = (rho.sum(axis = 1) / n)
       logrho_proj = np.log(rho_proj)
       hist = histogram.histogram(logrho_proj.ravel(), 0.01)
       write_histogram(hist, stats + "hist1D_logrho_y.txt")
       hist = histogram.histogram_weighted(logrho_proj.ravel(), rho_proj.ravel(), 0.01)
       write_histogram(hist, stats + "hist1D_logrho_y_weighted_M.txt")
       n = rho_proj.size
       s = logrho_proj.sum() / n
       s2 = (logrho_proj**2).sum() / n
       rho2D = rho_proj.sum() / n
       write_number(s, stats + "logrho_mean_2D_y.txt")
       write_number(np.sqrt(s2 - s**2), stats + "logrho_sigma_2D_y.txt")
       write_number(rho2D, stats + "rho_2D_y.txt")
       for power in range(1,max_power + 1):
              rho_pow = (rho_proj**power).sum() / n
              write_number(rho_pow, stats + "rho_pow_" + str(power) + "_2D_y.txt")
       
       #1D projections
       n = rho.shape[0]
       rho_proj1D = (rho_proj.sum(axis = 0) / n).ravel()
       n = rho_proj1D.size
       s = np.log(rho_proj1D).sum() / n
       s2 = (np.log(rho_proj1D)**2).sum() / n
       write_number(s, stats + "logrho_mean_1D_yx.txt")
       write_number(np.sqrt(s2 - s**2), stats + "logrho_sigma_1D_yx.txt")
       for power in range(1,max_power + 1):
              rho_pow = (rho_proj1D**power).sum() / n
              write_number(rho_pow, stats + "rho_pow_" + str(power) + "_1D_yx.txt")
       
       n = rho.shape[0]
       rho_proj1D = (rho_proj.sum(axis = 1) / n).ravel()
       n = rho_proj1D.size
       s = np.log(rho_proj1D).sum() / n
       s2 = (np.log(rho_proj1D)**2).sum() / n
       write_number(s, stats + "logrho_mean_1D_yz.txt")
       write_number(np.sqrt(s2 - s**2), stats + "logrho_sigma_1D_yz.txt")
       for power in range(1,max_power + 1):
              rho_pow = (rho_proj1D**power).sum() / n
              write_number(rho_pow, stats + "rho_pow_" + str(power) + "_1D_yz.txt")
       
       n = rho.shape[2]
       rho_proj = (rho.sum(axis = 2) / n)
       logrho_proj = np.log(rho_proj)
       hist = histogram.histogram(logrho_proj.ravel(), 0.01)
       write_histogram(hist, stats + "hist1D_logrho_z.txt")
       hist = histogram.histogram_weighted(logrho_proj.ravel(), rho_proj.ravel(), 0.01)
       write_histogram(hist, stats + "hist1D_logrho_z_weighted_M.txt")
       n = rho_proj.size
       s = logrho_proj.sum() / n
       s2 = (logrho_proj**2).sum() / n
       rho2D = rho_proj.sum() / n
       write_number(s, stats + "logrho_mean_2D_z.txt")
       write_number(np.sqrt(s2 - s**2), stats + "logrho_sigma_2D_z.txt")
       write_number(rho2D, stats + "rho_2D_z.txt")
       for power in range(1,max_power + 1):
              rho_pow = (rho_proj**power).sum() / n
              write_number(rho_pow, stats + "rho_pow_" + str(power) + "_2D_z.txt")
       
       #1D projections
       n = rho.shape[0]
       rho_proj1D = (rho_proj.sum(axis = 0) / n).ravel()
       n = rho_proj1D.size
       s = np.log(rho_proj1D).sum() / n
       s2 = (np.log(rho_proj1D)**2).sum() / n
       write_number(s, stats + "logrho_mean_1D_zx.txt")
       write_number(np.sqrt(s2 - s**2), stats + "logrho_sigma_1D_zx.txt")
       for power in range(1,max_power + 1):
              rho_pow = (rho_proj1D**power).sum() / n
              write_number(rho_pow, stats + "rho_pow_" + str(power) + "_1D_zx.txt")
       
       n = rho.shape[0]
       rho_proj1D = (rho_proj.sum(axis = 1) / n).ravel()
       n = rho_proj1D.size
       s = np.log(rho_proj1D).sum() / n
       s2 = (np.log(rho_proj1D)**2).sum() / n
       write_number(s, stats + "logrho_mean_1D_zy.txt")
       write_number(np.sqrt(s2 - s**2), stats + "logrho_sigma_1D_zy.txt")
       for power in range(1,max_power + 1):
              rho_pow = (rho_proj1D**power).sum() / n
              write_number(rho_pow, stats + "rho_pow_" + str(power) + "_1D_zy.txt")

       # bulk histograms
       rho = rho.ravel()
       logrho = np.log(rho)
       v2 = (vx**2 + vy**2 + vz**2).ravel()
       hist = histogram.histogram(logrho, 0.01)
       write_histogram(hist, stats + "hist1D_logrho.txt")
       hist = histogram.histogram_weighted(logrho, rho, 0.01)
       write_histogram(hist, stats + "hist1D_logrho_weighted_M.txt")
       hist = histogram.histogram_weighted(logrho, 0.5 * rho * v2, 0.01)
       write_histogram(hist, stats + "hist1D_logrho_weighted_KE.txt")
       n = rho.size
       s = logrho.sum() / n
       s2 = (logrho**2).sum() / n
       rho3D = rho.sum() / n
       mach = np.sqrt(v2.sum() / (3 * n))
       write_number(s, stats + "logrho_mean_3D.txt")
       write_number(np.sqrt(s2 - s**2), stats + "logrho_sigma_3D.txt")
       write_number(mach, stats + "mach_1D.txt")
       write_number(rho3D, stats + "rho_3D.txt")
       #for power in range(1,max_power + 1):
       #       rho_pow = (rho**power).sum() / n
       #       write_number(rho_pow, stats + "rho_pow_" + str(power) + "_3D.txt")

def moments_from_fourier2D(fourier, axis):
       Nbins = fourier.shape[0]
       Rbins = np.arange(Nbins)
       xx2,yy2=np.mgrid[0:Nbins,0:Nbins]
       xx3,yy3,zz3=np.mgrid[0:Nbins,0:Nbins,0:Nbins]
       r2D=np.sqrt(xx2**2+yy2**2)
       r3D=np.sqrt(xx3**2+yy3**2+zz3**2)
       digitized2D=np.digitize(r2D,Rbins)
       digitized3D=np.digitize(r3D,Rbins)
       fourierSq = fourier * np.conjugate(fourier)
       shellAvg = np.array([fourierSq[digitized2D == i].mean() if (digitized2D==i).any() else 0 for i in Rbins])
       radius2D = np.array([(digitized2D == i).sum() for i in Rbins])
       radius3D = np.array([(digitized3D == i).sum() for i in Rbins])

       rho1 = fourier[0,0]
       rho22D1 = (np.abs(fourier)**2).sum(axis = None)
       rho22D2 = (shellAvg * radius2D).sum(axis = None)
       rho23D = (shellAvg * radius3D).sum(axis = None)
       fourier_ext = fourier[:, :, np.newaxis, np.newaxis]
       i_ind, j_ind, k_ind, l_ind = np.indices((Nbins, Nbins, Nbins, Nbins))
       fourier_comb = fourier_ext[(i_ind + k_ind) % Nbins, (j_ind + l_ind) % Nbins]
       rho32D = (fourier_ext * fourier[np.newaxis, np.newaxis, :, :] * fourier_comb).sum()
       write_number(rho1, stats + "rho_pow_1_2D_fourier_" + axis + ".txt")
       write_number(rho22D1, stats + "rho_pow_2_2D_fourier_" + axis + ".txt")
       write_number(rho22D2, stats + "rho_pow_2_2D_fourier_" + axis + "_test.txt")
       write_number(rho23D, stats + "rho_pow_2_3D_fourier_" + axis + ".txt")
       write_number(rho32D, stats + "rho_pow_3_2D_fourier_" + axis + ".txt")

def fourier2D(rho, axis = "x"):
       n = rho.shape[0]
       ax = 0
       if (axis == "y"):
              n = rho.shape[1]
              ax = 1
       elif (axis == "z"):
              n = rho.shape[2]
              ax = 2
       print("Computing Fourier image of density projection along " + axis)
       rho_proj = (rho.sum(axis = ax) / n)
       n = rho_proj.size
       fourier = np.fft.fft2(rho_proj) / n
       re = np.real(fourier)
       im = np.imag(fourier)
       kx = 2 * np.pi * np.fft.fftfreq(rho_proj.shape[0])
       ky = 2 * np.pi * np.fft.fftfreq(rho_proj.shape[1])
       write_fourier2D(kx, ky, re, stats + "rho_re_fourier_2D_" + axis + ".txt")
       write_fourier2D(kx, ky, im, stats + "rho_im_fourier_2D_" + axis + ".txt")
       moments_from_fourier2D(fourier, axis)

mach =   "1"
alfven = "1"

sims = ["1_1", "1_2", "1_half", "2_1", "2_2", "2_half", "3_1", "3_2", "3_half", "4_1", "4_2", "4_half", "5_1", "5_2", "5_half", "6_1", "6_2", "6_half"]
frames = [31, 31, 31, 31, 31, 85, 76, 40, 92, 51, 51, 44, 48, 59, 59, 51, 51, 45]

for sim_index in range(0, len(sims)):
       sim = sims[sim_index]
       frame = frames[sim_index]
       print("Simulation: %s"%sim)

       pics = "./results/" + sim + "/pics/"
       stats = "./results/" + sim + "/stats/"

       if not os.path.exists("./results/" + sim):
              os.mkdir("./results/" + sim)

       if not os.path.exists("./results/" + sim + "/pics/"):
              os.mkdir("./results/" + sim + "/pics/")

       if not os.path.exists("./results/" + sim + "/stats/"):
              os.mkdir("./results/" + sim + "/stats/")

       ds = yt.load("/data/cb1/Projects/P49_EE_BB/%s/DD%04d/data%04d"%(sim, frame, frame))
       #ds = yt.load("/data/cb1/Projects/P49_EE_BB/%s/DD0011/data0011"%(sim))
       cg = ds.covering_grid(0, [0.0]*3, [512]*3)
       rho = cg["density"].v
       #fourier2D(rho, axis="x")
       #fourier2D(rho, axis="y")
       #fourier2D(rho, axis="z")
       rhoX = hessian.extend_periodic(rho.sum(axis=0),1) / 512
       rhoY = hessian.extend_periodic(rho.sum(axis=1),1) / 512
       rhoZ = hessian.extend_periodic(rho.sum(axis=2),1) / 512
       logrho_periodic = hessian.extend_periodic(np.log(rho),1)
       Bx = cg["magnetic_field_x"].v
       By = cg["magnetic_field_y"].v
       Bz = cg["magnetic_field_z"].v
       
       #vx = cg["x-velocity"].v
       #vy = cg["y-velocity"].v
       #vz = cg["z-velocity"].v

       print("Data loaded!\n")
       #density_stats(rho, vx, vy, vz)
       
       dds = cg.dds

       '''
       stats_Bv(Bx, By, Bz, vx, vy, vz, rho)
       '''
       t = time.perf_counter()
       
       hess = hessian.hessian(logrho_periodic,dds)
       print("time to obtain hessian: %ds"%(time.perf_counter() - t))
       t = time.perf_counter()

       e,v=hessian.eigensystem(hess)
       print("time to obtain eigensystem: %ds"%(time.perf_counter() - t))
       t = time.perf_counter()

       hessian.sort_hessian(e,v)
       print("time to sort the eigensystem: %ds"%(time.perf_counter() - t))

       t = time.perf_counter()

       '''
       print("Extracting 3D statistics...")
       blob_statistics(e,rho)
       filament_statistics_3D(e,v, rho, Bx, By, Bz)
       #sheet_statistics(e,v,rho)
       '''
       '''
       print("Extracting 2D statistics...")
       hess = hessian.hessian(rhoX,dds[[1,2]])
       e,v=hessian.eigensystem(hess)
       hessian.sort_hessian(e,v)
       filament_statistics_2D(e,v, rho, Bx, By, Bz, "x")

       hess = hessian.hessian(rhoY,dds[[0,2]])
       e,v=hessian.eigensystem(hess)
       hessian.sort_hessian(e,v)
       filament_statistics_2D(e,v, rho, Bx, By, Bz, "y")

       hess = hessian.hessian(rhoZ,dds[[0,1]])
       e,v=hessian.eigensystem(hess)
       hessian.sort_hessian(e,v)
       filament_statistics_2D(e,v, rho, Bx, By, Bz, "z")
       '''

       #mask = mask3D_smooth(e, 1)
       n = rho.size
       eps = 0.01
       histE = histogram.histogram2D((e[...,0]).ravel(), (e[...,1]).ravel(), 1000, 1000)
       write_histogram(histE, stats + "eigenvalue_12.txt")
       histE = histogram.histogram2D((e[...,0]).ravel(), (e[...,2]).ravel(), 1000, 1000)
       write_histogram(histE, stats + "eigenvalue_13.txt")
       histE = histogram.histogram2D((e[...,1]).ravel(), (e[...,2]).ravel(), 1000, 1000)
       write_histogram(histE, stats + "eigenvalue_23.txt")
       '''
       meanE = (e[...,0]).sum() / n
       sigmaE = np.sqrt((e[...,0]**2).sum() / n - meanE)
       write_number(meanE, stats + "eigenvalue_1_mean.txt")
       write_number(sigmaE, stats + "eigenvalue_1_sigma.txt")
       '''
       meanE = (e[...,1]).sum() / n
       sigmaE = np.sqrt((e[...,1]**2).sum() / n - meanE)
       #write_number(meanE, stats + "eigenvalue_2_mean.txt")
       #write_number(sigmaE, stats + "eigenvalue_2_sigma.txt")
       '''
       meanE = (e[...,2]).sum() / n
       sigmaE = np.sqrt((e[...,2]**2).sum() / n - meanE)
       write_number(meanE, stats + "eigenvalue_3_mean.txt")
       write_number(sigmaE, stats + "eigenvalue_3_sigma.txt")
       '''
       mask = mask3D_smooth3(e, 2 * sigmaE)
       postfix = "_smooth_alt_2"
       
       for axis in ["x", "y", "z"]:
              psi_B, q, u = getQU(rho, Bx, By, Bz, ax=axis)
              kx, ky, fourierE, fourierB, modeE, modeB = getEB(q, u, axis, "")
              write_fourier2D(kx, ky, np.real(fourierE), stats + "mode_re_E_fourier_" + axis + ".txt")
              write_fourier2D(kx, ky, np.imag(fourierE), stats + "mode_im_E_fourier_" + axis + ".txt")
              write_fourier2D(kx, ky, np.real(fourierB), stats + "mode_re_B_fourier_" + axis + ".txt")
              write_fourier2D(kx, ky, np.imag(fourierB), stats + "mode_im_B_fourier_" + axis + ".txt")

              write_fourier2D(kx, ky, np.real(modeE), stats + "mode_re_E_" + axis + ".txt")
              write_fourier2D(kx, ky, np.imag(modeE), stats + "mode_im_E_" + axis + ".txt")
              write_fourier2D(kx, ky, np.real(modeB), stats + "mode_re_B_" + axis + ".txt")
              write_fourier2D(kx, ky, np.imag(modeB), stats + "mode_im_B_" + axis + ".txt")

              psi_B, q, u = getQU(rho * mask, Bx, By, Bz, ax=axis)
              kx, ky, fourierE, fourierB, modeE, modeB = getEB(q, u, axis, "_masked" + postfix)
              write_fourier2D(kx, ky, np.real(fourierE), stats + "mode_re_E_fourier_" + axis + "_masked" + postfix + ".txt")
              write_fourier2D(kx, ky, np.imag(fourierE), stats + "mode_im_E_fourier_" + axis + "_masked" + postfix + ".txt")
              write_fourier2D(kx, ky, np.real(fourierB), stats + "mode_re_B_fourier_" + axis + "_masked" + postfix + ".txt")
              write_fourier2D(kx, ky, np.imag(fourierB), stats + "mode_im_B_fourier_" + axis + "_masked" + postfix + ".txt")

              write_fourier2D(kx, ky, np.real(modeE), stats + "mode_re_E_" + axis + "_masked" + postfix + ".txt")
              write_fourier2D(kx, ky, np.imag(modeE), stats + "mode_im_E_" + axis + "_masked" + postfix + ".txt")
              write_fourier2D(kx, ky, np.real(modeB), stats + "mode_re_B_" + axis + "_masked" + postfix + ".txt")
              write_fourier2D(kx, ky, np.imag(modeB), stats + "mode_im_B_" + axis + "_masked" + postfix + ".txt")

              psi_B, q, u = getQU(rho * (1 - mask), Bx, By, Bz, ax=axis)
              kx, ky, fourierE, fourierB, modeE, modeB = getEB(q, u, axis, "_complement" + postfix)
              write_fourier2D(kx, ky, np.real(fourierE), stats + "mode_re_E_fourier_" + axis + "_complement" + postfix + ".txt")
              write_fourier2D(kx, ky, np.imag(fourierE), stats + "mode_im_E_fourier_" + axis + "_complement" + postfix + ".txt")
              write_fourier2D(kx, ky, np.real(fourierB), stats + "mode_re_B_fourier_" + axis + "_complement" + postfix + ".txt")
              write_fourier2D(kx, ky, np.imag(fourierB), stats + "mode_im_B_fourier_" + axis + "_complement" + postfix + ".txt")

              write_fourier2D(kx, ky, np.real(modeE), stats + "mode_re_E_" + axis + "_complement" + postfix + ".txt")
              write_fourier2D(kx, ky, np.imag(modeE), stats + "mode_im_E_" + axis + "_complement" + postfix + ".txt")
              write_fourier2D(kx, ky, np.real(modeB), stats + "mode_re_B_" + axis + "_complement" + postfix + ".txt")
              write_fourier2D(kx, ky, np.imag(modeB), stats + "mode_im_B_" + axis + "_complement" + postfix + ".txt")