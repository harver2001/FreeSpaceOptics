import numpy as np
from scipy.special import erf,gamma
from scipy.constants import Boltzmann
from scipy.constants import Planck
import math
import matplotlib.pyplot as plt

N = 35
Hm = 3.2
var_to = 1e-3
var_ro = 1e-3
var_tp = 30e-2
var_rp = 30e-2
var_R = 2
ra = 4e-2
wz = 1
hl_dB = 10
hl = 10**(-hl_dB/10)
e =1.6e-19
G = 15
Be = 1e9
Nb = 1e-3
Lamda = 1550e-9
Fov = 5e-3 
F = 2.75
K= 1.38e10-23
Pb = np.pi**2 *ra*Be*Nb*Lamda*Fov /4
hl = 10**(-hl_dB/10)
# Z = 500
Z = np.array([300,500,700])
FOV = np.arange(0.01,0.1,0.01)
v = np.sqrt(np.pi/2)*(ra/wz)
wz_eq = wz**2*(np.sqrt(np.pi)*erf(v))/(2*v*np.exp(-v**2))
A0 = (erf(v))**2
T = (wz_eq)/(4*((Z*var_to)**2 + var_tp**2 + var_rp**2)) + 1
C1 = (T/((A0*hl)*T)) * np.exp(0.5*var_R*T*(1+T))
C2 = T - 1
C3 = 0.5 * var_R * (1 + 2*T)
Pt= 316.2
Rl = 1e3 # Load Resistance
Tr = 300 # Receiver Temperature
var_th = ((4*Boltzmann)*Tr*(1e9)) / (Rl)
var_bg = 2*e*G*F*Be*Nb*Pb
var_n = var_bg + var_th
Mew = 0.8
c=3e8
v1= c/(Lamda)
eta = e*G*Mew/(Planck*v1)
b1 = eta*Pt/(np.sqrt(8*(var_th+var_bg)))
Hm = 3.2
N = 35
YTh = 1e6
R=1

Pt = np.array([10,20,30,40,50])
Hth=[0] * len(Pt)

for i in range(len(Pt)):
    Hth[i] = np.sqrt(var_n * YTh) / (R * Pt[i])

# alpha = 4
# beta = 1.9
alpha = 1/((np.exp(0.49*var_R )/(1+1.11* (var_R)**(12/5))**(7/6)))
beta =  1/((np.exp(0.51*var_R )/(1+0.69* (var_R)**(12/5))**(5/6)))

Tow = wz_eq**2 /(4*(Z**2 *(var_to**2)+ var_tp**2+ var_rp**2))

P1 = 0

Pout = [[0] * len(Pt) for _ in range(len(Z))]
for x in range(len(Z)):
    for j in range (len(Pt)):
        for i in range (0,N+1):
            c4 = (np.pi*T[x])/(((hl*A0)**T[x])*gamma(alpha)*gamma(beta)*np.sin(np.pi*(alpha - beta)))
            c5 = ((alpha*beta)**(i + beta))/((i + beta - T[x])*gamma(i - alpha + beta + 1)*math.factorial(i))
            c6 = ((alpha*beta)**(i + alpha))/((i + alpha - T[x])*gamma(i - alpha + beta + 1)*math.factorial(i))
            c1 = (c4*c5*Hm**(i + beta - T[x]) - c4*c6*Hm**(i + alpha - T[x]))
            c2 = (c4*c6)/((hl*A0)**(i + alpha - T[x]))
            c3 = (c4*c5)/((hl*A0)**(i + beta - T[x]))
            P1 = P1+(c1/Tow[x])*((Hth[j])**(Tow[x])) + ((c2/(N+ alpha))*((Hth[j])**(N+alpha))) - ((c3/(N+beta))*((Hth[j])**(N+beta)))
        Pout[x][j] = b1+(1-b1)*P1

# Normalize Pout values
Pout_normalized = [[Pout[i][j] / np.linalg.norm(Pout[i]) for j in range(len(Pout[i]))] for i in range(len(Pout))]

# Plot all normalized Pout values in a single graph
plt.figure(figsize=(10, 6))
for x in range(len(Z)):
    plt.semilogy(Pt, Pout_normalized[x], label=f'Z={Z[x]}')

plt.xlabel('Transmitted Power (dBm)')
plt.ylabel('Outage Probability (Pout)')
plt.legend()
plt.grid(True)
plt.show()
