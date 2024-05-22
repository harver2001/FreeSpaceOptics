import numpy as np
from scipy.special import erf, gamma
from scipy.constants import Boltzmann, Planck
import math
import matplotlib.pyplot as plt

# Constants
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
e = 1.6e-19
G = 15
Be = 1e9
Nb = 1e-3
Lamda = 1550e-9
F = 2.75
K = 1.38e-23
Z = 500
Rl = 1e3  # Load Resistance
Tr = 300  # Receiver Temperature
var_th = ((4*Boltzmann)*Tr*(1e9)) / (Rl)
Mew = 0.8
c = 3e8
v1 = c / (Lamda)
eta = e*G*Mew / (Planck*v1)
YTh = 1e6
R = 1
Pt = np.array([0.01, 0.025 , 0.045, 0.06])

FOVs = [5e-3, 5.5e-3, 6e-3, 6.5e-3]

for Fov in FOVs:
    turn = 1
    Pb = np.pi**2 * ra * Be * Nb * Lamda * Fov / 4
    var_bg = 2 * e * G * F * Be * Nb * Pb
    var_n = var_bg + var_th
    b1 = eta * Pt / (np.sqrt(8 * (var_th + var_bg)))
    Hth = [np.sqrt(var_n * YTh) / (R * p) for p in Pt]

    alpha = 1 / ((np.exp(0.49 * var_R) / (1 + 1.11 * (var_R)**(12/5))**(7/6)))
    beta = 1 / ((np.exp(0.51 * var_R) / (1 + 0.69 * (var_R)**(12/5))**(5/6)))

    v = np.sqrt(np.pi / 2) * (ra / wz)
    wz_eq = wz**2 * (np.sqrt(np.pi) * erf(v)) / (2 * v * np.exp(-v**2))
    Tow = wz_eq**2 / (4 * (Z**2 * (var_to**2) + var_tp**2 + var_rp**2))
    A0 = (erf(v))**2
    T = (wz_eq) / (4 * ((Z * var_to)**2 + var_tp**2 + var_rp**2)) + 1

    Pout = []
    for j in range(len(Pt)):
        P1 = 0
        for i in range(N + 1):
            c4 = (np.pi * T) / (((hl * A0)**T) * gamma(alpha) * gamma(beta) * np.sin(np.pi * (alpha - beta)))
            c5 = ((alpha * beta)**(i + beta)) / ((i + beta - T) * gamma(i - alpha + beta + 1) * math.factorial(i))
            c6 = ((alpha * beta)**(i + alpha)) / ((i + alpha - T) * gamma(i - alpha + beta + 1) * math.factorial(i))
            c1 = (c4 * c5 * Hm**(i + beta - T)) - (c4 * c6 * Hm**(i + alpha - T))
            P1 += (c1 / Tow) * ((Hth[j])**Tow)
        Pout.append((b1 + (1 - b1) * P1)*0.1)
    Pout = Pout / np.max(Pout)
    plt.plot(Pt * 1000, Pout*0.1)

plt.xlabel('Power Transmitted (dBm)')
plt.ylabel('Outage Probability')
# plt.title('Outage Probability vs Power Transmitted for Different FOVs')
# plt.ylim([10**-5, 1.1])
plt.legend()
plt.show()