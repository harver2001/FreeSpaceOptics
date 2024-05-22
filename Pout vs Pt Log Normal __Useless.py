import numpy as np
from scipy.special import erf
from scipy.constants import Boltzmann
from scipy.constants import Planck
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
Z = 500
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

Pt = np.array([10,20,30,40,50,60,70])
Hth=[0] * len(Pt)

for i in range(len(Pt)):
    Hth[i] = np.sqrt(var_n * YTh) / (R * Pt[i])

alpha = 1/((np.exp(0.49*var_R )/(1+1.11* (var_R)**(12/5))**(7/6))-1)
beta =  1/((np.exp(0.51*var_R )/(1+0.69* (var_R)**(12/5))**(5/6))-1)

Tow = wz_eq**2 /(4*(Z**2 *(var_to**2)+ var_tp**2+ var_rp**2))

Pout = [0] * len(Pt)

for j in range(len(Pt)):
    P1 = 0 
    for i in range(0, N+1):
        P1 += (C1/Tow)*((Hth[j])**(Tow)) + ((C2/(N+ alpha))*((Hth[j])**(N+alpha))) - ((C3/(N+beta))*((Hth[j])**(N+beta)))
    Pout[j] = b1+(1-b1)*P1

plt.plot(Pt, Pout)
plt.xlabel('Pt')
plt.ylabel('Poutage')
plt.show()