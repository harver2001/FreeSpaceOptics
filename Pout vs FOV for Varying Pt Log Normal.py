import numpy as np
from scipy.special import erf
from scipy.stats import norm
from scipy.constants import Boltzmann, Planck
import matplotlib.pyplot as plt
import math

q_charge = 1.602e-19  # C
gain_apd_values = np.array([15])
f_apd = 2.75

h = np.linspace(0,1e-3,50)
ra = 4e-2
wz = 1
var_R = 0.2
var_to = 1e-3
var_ro = 1e-3
var_tp = 30e-2
var_rp = 30e-2
hl_dB = 10
hl = 10**(-hl_dB/10)
Z = 500
Z1 = [500, 1000, 1500, 2000, 2500, 3000]
FOV = np.linspace(0.005,0.05,100)
Theta_FOV = 5e-3 
v = np.sqrt(np.pi/2)*(ra/wz)
wz_eq = wz**2*(np.sqrt(np.pi)*erf(v))/(2*v*np.exp(-v**2))
A0 = (erf(v))**2
T = (wz_eq)/(4*((Z*var_to)**2 + var_tp**2 + var_rp**2)) + 1
C1 = (T/((A0*hl)*T)) * np.exp(0.5*var_R*T*(1+T))
C2 = T - 1
C3 = 0.5 * var_R * (1 + 2*T)
a=[5/24,4/24,1/24]
az = [2,11/20,1/2]
Rl = 1e3 # Load Resistance
Tr = 300 # Receiver Temperature
lamda = 1550e-9
c = 3e8
v1 = c/lamda
G = 30 # Average APD Gain
Mew = 0.8
Be = 1e9 # electrical bandwidth
B0 = 1e-9 # optical filter bandwidth
Nb = 1e-3 # spectral radiance of background radiation 
eta = q_charge*G*Mew/(Planck*v1)

Pb = (np.pi**2 * ra**2 * B0 * Nb * FOV**2) / 4
var_th = (4*Boltzmann*Tr*(1e9)) / (Rl)
var_bg = 2*q_charge*G*f_apd*Be*eta*Pb
var_n = var_bg + var_th
C4 = []

for i in range(len(az)):
    x = -az[i]*np.log(A0*hl)+az[i]*C3 - var_R*(C2+1)/2
    C4.append(x)
    
i=0
j=0
Pt = np.arange(1,30,5)
YTh = 1e6
R = 1
Hth = np.zeros(len(Pt))

for i in range(len(Pt)):
    Hth[i] = (i+1.5)*1e-7


plt.figure()
for h in range(len(Hth)):
    Pout = []
    for i in range(len(FOV)):
        temp = 0
        for j in range(len(a)):
            temp += (a[j]/2)*np.sqrt(np.pi/az[j])*np.exp((C4[j]**2 - (az[j]*(np.log(A0*hl) - C3))**2)/az[j])*(1 + erf(np.sqrt(az[j])*np.log(Hth[h]) + C4[j]/np.sqrt(az[j])))
        temp = (Hth[h]**(C2 + 1))/(C2 + 1) - temp
        temp = temp*C1*(1 - np.exp(-(FOV[i]**2)/(2*(var_to**2 + var_ro**2)))) + np.exp(-(FOV[i]**2)/(2*(var_to**2 + var_ro**2)))
        Pout.append(temp)
    plt.semilogy(FOV*1000, Pout, label=f'Pt={Pt[len(Pt)-1-h]} dBm')
plt.xlabel('FOV(mrad)')
plt.ylabel('Outage Probablity')
plt.xlim([7,12])
plt.ylim([10**-7.5,10**-5])
plt.legend()
plt.show()

