import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.constants import Boltzmann, Planck
from scipy.special import erf

# Constants and Parameters Initialization
q_charge = 1.602e-19  # Charge of an electron in Coulombs
gain_apd_values = np.array([15])
f_apd = 2.75  # Excess noise factor of APD
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
FOV = 5e-3 
lamda = 1550e-9
c = 3e8
v = np.sqrt(np.pi/2)*(ra/wz)
v1 = c/lamda
wz_eq = wz**2*(np.sqrt(np.pi)*erf(v))/(2*v*np.exp(-v**2))
A0 = (erf(v))**2
T = (wz_eq)/(4*((Z*var_to)**2 + var_tp**2 + var_rp**2)) + 1
C1 = (T/((A0*hl)*T)) * np.exp(0.5*var_R*T*(1+T))
C2 = T - 1
C3 = 0.5 * var_R * (1 + 2*T)
Mew = 0.8
Be = 1e9  # Electrical bandwidth
B0 = 1e-9  # Optical filter bandwidth
Nb = 1e-3  # Spectral radiance of background radiation 
G = 30  # Average APD Gain
eta = q_charge*G*Mew/(Planck*v1)
Pb = (np.pi**2 * ra**2 * B0 * Nb * FOV**2) / 4
Rl = 1e3  # Load Resistance
Tr = 300  # Receiver Temperature
var_th = (4*Boltzmann*Tr*(1e9)) / (Rl)
var_bg = 2*q_charge*G*f_apd*Be*eta*Pb
a = [5/24, 4/24, 1/24]
az = [2, 11/20, 1/2]
C4 = []

for i in range(len(az)):
    x = -az[i]*np.log(A0*hl)+az[i]*C3 - var_R*(C2+1)/2
    C4.append(x)

Pt = np.linspace(-0.010, 0.050, 1000) / 1000 
b1 = (eta * Pt) / np.sqrt(8 * (var_th + var_bg))
h1 = A0*hl*np.exp(-1*C3)
Pe = np.zeros(len(Pt))

for i in range(len(Pt)):
    Pe[i] = 0.5 * np.exp((-1 * FOV**2) / (2 * (var_to**2 + var_ro**2)))
    temp = 0
    for j in range(3):
        temp += (0.5 * a[j] * h1**(C1 + 1)) * ((1 / (C2 + 1)) + np.sqrt((var_R * np.pi) / az[j]) * (np.exp(var_R * (C2 + 1)**2) / (4 * az[j])) * erf(((np.sqrt(var_R)) * (C2 + 1)) / (2 * np.sqrt(az[j]))))
        val1 = 2 / np.sqrt(np.pi)
        val = 0
        for m in range(10):
            val += ((-1)**m * (b1[i] * h1)**(2 * m + 1)) / (math.factorial(m) * (2 * m + 1)) * ((-1 / (C2 + 2 * m + 2)) + np.sqrt((var_R * np.pi) / az[j]) * np.exp((var_R * (C2 + 2 * m + 2)**2) / (4 * az[j])) * erf((np.sqrt(var_R) * (C2 + 2 * m + 2)) / (2 * np.sqrt(az[j]))))
        val *= val1
        temp += val
    Pe[i] += C1 * (1 - np.exp((-1 * FOV**2) / (2 * (var_to**2 + var_ro**2)))) * temp

Pe_scaled = 10**(-10) + (1 - 10**(-10)) * (Pe - np.min(Pe)) / ((np.max(Pe) - np.min(Pe)) * 1000)

# Plotting
plt.figure()
plt.plot(Pt*1000, Pe_scaled , '-')  
plt.xlabel('Transmitted Power')
plt.ylabel('Bit Error Rate (BER)')
plt.xlim([0.010, 0.05])
plt.title('BER vs. Transmitted Power')
plt.grid(True)
plt.show()