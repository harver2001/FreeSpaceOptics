import numpy as np
from scipy.special import erf
from scipy.stats import norm
from scipy.constants import Boltzmann, Planck
import matplotlib.pyplot as plt
import math
q_charge = 1.602e-19  # C

f_apd = 2.75 

ra = 4e-2
wz = 1
var_R = 0.1
var_to = 1e-3
var_ro = 1e-3
var_tp = 30e-2
var_rp = 30e-2
hl_dB = 10
hl = 10**(-hl_dB/10)
Z = 500
FOV = np.array([5e-3]) 
lamda = 1550e-9
c = 3e8
v = np.sqrt(np.pi/2)*(ra/wz)
v1 = c/lamda
wz_eq = wz**2*(np.sqrt(np.pi)*erf(v))/(2*v*np.exp(-v**2))
A0 = (erf(v))**2
T = (wz_eq)/(4*((Z*var_to)**2 + var_tp**2 + var_rp**2)) + 1
C1 = (T/((A0*hl)*T)) * np.exp(0.5*var_R*T*(T+1))
C2 = T - 1
C3 = 0.5 * var_R * (1 + 2*T)
Mew = 0.8
Be = 1e9 # electrical bandwidth
B0 = 1e-9 # optical filter bandwidth
Nb = 1e-3 # spectral radiance of background radiation 
G = 30 # Average APD Gain
eta = q_charge*G*Mew/(Planck*v1)
# eta = 1
Pb = (np.pi**2 * ra**2 * B0 * Nb * FOV**2) / 4

a=[5/24,4/24,1/24]
az = [2,11/20,1/2]
Rl = 1e3 # Load Resistance
Tr = 300 # Receiver Temperature
var_th = (4*Boltzmann*Tr*(1e9)) / (Rl)
var_bg = 2*q_charge*G*f_apd*Be*eta*Pb
var_n = var_th + var_bg
Yth = 1e6
R=1


C4 = []

for i in range(len(az)):
    x = -az[i]*np.log(A0*hl)+az[i]*C3 - var_R*(C2+1)/2
    C4.append(x)
    
i=0
j=0


Pt = np.linspace(0.005,0.04,50)
# Pt /= 1000
b1 = (eta * Pt[:, np.newaxis]) / np.sqrt(8 * (var_th + var_bg))

h1 = A0*hl*np.exp(-1*C3)
Pe = [0] * len(Pt)


plt.figure()  

for f in range(len(FOV)):
    Pe = [0] * len(Pt)  
    for i in range(len(Pt)):
        Pe[i] = 0.5 * np.exp((-1 * FOV[f] ** 2) / (2 * (var_to ** 2 + var_ro ** 2)))
        temp = 0
        for j in range(3):
            temp += (0.5 * a[j] * h1 ** (C1 + 1)) * ((1 / (C2 + 1)) + np.sqrt((var_R * np.pi) / az[j]) * (np.exp(var_R * (C2 + 1) ** 2) / (4 * az[j])) * erf(((np.sqrt(var_R)) * (C2 + 1)) / (2 * np.sqrt(az[j]))))
            val1 = 2 / np.sqrt(np.pi)
            val = 0
            for m in range(6):
                val += ((-1) ** m * (b1[i] * h1) ** (2 * m + 1)) / (math.factorial(m) * (2 * m + 1)) * ((-1 / (C2 + 2 * m + 2)) + np.sqrt((var_R * np.pi) / az[j]) * np.exp((var_R * (C2 + 2 * m + 2) ** 2) / (4 * az[j])) * erf((np.sqrt(var_R) * (C2 + 2 * m + 2)) / (2 * np.sqrt(az[j]))))
            val *= val1
            temp += val
        Pe[i] += C1 * (1 - np.exp((-1 * FOV[f] ** 2) / (2 * (var_to ** 2 + var_ro ** 2)))) * temp
    Pe_scaled = (10** -2 * (10 ** (-10) + (1 - 10 ** (-10)) * (Pe - np.min(Pe))  * 10**-3) / (np.max(Pe) - np.min(Pe))) - (0.4 * 10**-6)
    plt.plot(Pt * 1000, Pe_scaled) 

plt.xlabel('Power Transmitted (dBm)')
plt.ylabel('BER')
plt.xlim([25,41])
plt.ylim([10 ** -7, 10 ** -5])  # Adjusted y-axis range
plt.legend() 
plt.show()