import math
import numpy as np
from scipy.special import erf, gamma
from scipy.constants import Boltzmann, Planck, speed_of_light
import matplotlib.pyplot as plt

# Constants and parameters
Pt_dB = np.linspace(-10, 30, 100)
Pt_lin = 10**(Pt_dB/10)
wavelength = 1550e-9
gain_APD = 15
mu_mod = 0.8
Tr = 300
Be = 1e9
Bo = 1e-9
Nb_lambda = 1e7
FOV = np.array([12e-3, 50e-3])
F_APD = 2.75

var_R = 2
hm = 3.5
N = 10

R_load = 1e3
hl_dB = -20
hl_lin = 10**(hl_dB/10)
ra = 4e-2
wz = 1.35  # Important
Z = 500
sd_to = 1e-3
sd_ro = 1e-3
sd_tp = 30e-2
sd_rp = 30e-2

v = np.sqrt(np.pi/2)*(ra/wz)
A0 = (erf(v))**2
wz_eq_sq = (wz**2)*(np.sqrt(np.pi/4)*erf(v))/(v*np.exp(-v**2))
T = wz_eq_sq/(4*((Z*sd_to)**2 + sd_tp**2 + sd_rp**2))

C1 = (T/((A0*hl_lin)**T))*np.exp(0.5*T*(1 + T)*var_R**2)
C3 = 0.5*var_R*(1 + 2*T)
h1 = A0*hl_lin*np.exp(-C3)

alpha = (np.exp((0.49*var_R)/((1+1.11*var_R**(6/5))**(7/6))) - 1)**(-1)
beta = (np.exp((0.51*var_R)/((1+0.69*var_R**(6/5))**(7/6))) - 1)**(-1)

n_eta = (1.6e-19)*gain_APD*mu_mod/((6.62e-34)*speed_of_light/wavelength)
Pb = ((np.pi*ra/2)**2)*Bo*Nb_lambda*(FOV**2)
var_th = 4*Boltzmann*Tr*Be/R_load
var_b = 2*(1.6e-19)*gain_APD*F_APD*Be*n_eta*Pb
b1 = np.zeros((len(FOV), len(Pt_lin)))

for i in range(len(FOV)):
    for j in range(len(Pt_lin)):
        b1[i, j] = (n_eta*Pt_lin[j])/np.sqrt(8*(var_th + var_b[i]))

BER = np.zeros((len(FOV), len(b1[0])))

for i in range(len(FOV)):
    for j in range(len(Pt_lin)):
        for n in range(N):
            c4 = (np.pi*T)/(((hl_lin*A0)**T)*gamma(alpha)*gamma(beta)*np.sin(np.pi*(alpha - beta)))
            c5 = ((alpha*beta)**(n + beta))/((n + beta - T)*gamma(n - alpha + beta + 1)*math.factorial(n))
            c6 = ((alpha*beta)**(n + alpha))/((n + alpha - T)*gamma(n - alpha + beta + 1)*math.factorial(n))
            c1 = (c4*c5*hm**(n + beta - T) - c4*c6*hm**(n + alpha - T))
            c2 = (c4*c6)/((hl_lin*A0)**(n + alpha - T))
            c3 = (c4*c5)/((hl_lin*A0)**(n + beta - T))
            for m in range(1, 11):
                BER[i, j] += (((-1)**m)*(b1[i, j]**(2*m + 1))/((2*m + 1)*math.factorial(m)))*((c1*(hl_lin*A0*hm)**(T + 2*m + 1))/(T + 2*m + 1) + (c2*(hl_lin*A0*hm)**(n + alpha + 2*m + 1))/(n + alpha + 2*m + 1) - (c3*(hl_lin*A0*hm)**(n + beta + 2*m + 1))/(n + beta + 2*m + 1))
            BER[i, j] = -np.sqrt(4/np.pi)*BER[i, j]
            BER[i, j] += ((c1*(hl_lin*A0*hm)**T)/T) + ((c2*(hl_lin*A0*hm)**(n+alpha))/(n+alpha)) - ((c3*(hl_lin*A0*hm)**(n+beta))/(n+beta))
        BER[i, j] = BER[i, j]*C1*(1 - np.exp(-(FOV[i])**2/(2*(sd_to**2 + sd_ro**2))))
        BER[i, j] = BER[i, j] + np.exp(-(FOV[i])**2/(2*(sd_to**2 + sd_ro**2)))*(1/4)

print(BER[0])
plt.plot(Pt_dB, BER[0])
plt.show()
