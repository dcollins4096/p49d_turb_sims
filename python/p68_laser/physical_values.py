from dtools.starter1 import *
import unyt
unit_kappa = unyt.cm**2/unyt.g
unit_density = unyt.g/unyt.cm**3
alpha_B = 0.777*unit_kappa #cm^2/g
alpha_U = 4.697*unit_kappa 
alpha_D = 46.087*unit_kappa 
L_B = 0.06*unyt.cm #cm, 2 layers of 0.03
L_U = 0.14*unyt.cm
L_D = 0.06*unyt.cm
L_F = L_U+L_D
rho_B = 1.845*(unit_density) #g/cm^3
rho_U = 0.090*(unit_density)
rho_D = 0.099*(unit_density)
r = rho_D/rho_U
tau_B = rho_B*L_B*alpha_B
tau_U = rho_U*L_U*alpha_U
tau_D = rho_D*L_D*alpha_D
t1 = 1/(alpha_D*L_D*r+alpha_U*L_U)
t2 = np.exp(alpha_B*L_B*rho_B)

def compute_I0(I_pixel):
    I_0 = I_pixel * np.exp(tau_B+tau_U+tau_D)
    return I_0
def compute_rho(I_pixel, I_0):
    rho =  (np.log(I_0/I_pixel)-tau_B)*t1
    return rho




