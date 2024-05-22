import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
from scipy.constants import Boltzmann

# Constants
q_charge = 1.602e-19  # C
m_mod = 1
R_load = 1e3  # Ohm
Fn = 2
Rb = 2e9
kA = 0.7
R = 1
D = 0.02  # m
beta_v = (0.1 * np.log(10)) / 10000  # dB/km ---> 1/km
phi_divergence = 1e-3  # rad
wav_lambda = 1.55e-6  # m
T_Rx = 300  # K
Cn2 = 6e-14
Ps = np.logspace(-1, 7, 20)  # W
Ps_dB = 10 * np.log10(Ps)
L = np.array([1e3, 1.3e3, 1.6e3])  # m
w20 = np.array([2.22939364554e-13, 4.39934099226e-10, 1.08606937077e-7, 7.8025564785e-6, 0.000228338636017, 0.00324377334224, 0.0248105208875, 0.10901720602, 0.286675505363, 0.462243669601, 0.462243669601, 0.286675505363, 0.10901720602, 0.0248105208875, 0.00324377334224, 0.000228338636017, 7.8025564785e-6, 1.08606937077e-7, 4.39934099226e-10, 2.22939364554e-13])
x20 = np.linspace(-5.38748089001, 5.38748089001, 20)
gain_apd = 15
delta_f = Rb / 2
a = ((D / (phi_divergence * L))**2) * np.exp(-beta_v * L)

FA = kA * gain_apd + (1 - kA) * (2 - (1 / gain_apd))
var_R = 1.23 * Cn2 * ((2 * np.pi / wav_lambda)**(7/6)) * (L**(11/6))
d = np.sqrt(((2 * np.pi / wav_lambda) * D**2) / (4 * L))
var_s = np.exp((0.49 * var_R) / ((1 + 0.18 * d**2 + 0.56 * var_R**(6/5))**(7/6)) + (0.51 * var_R) / ((1 + 0.9 * d**2 + 0.62 * (d**2) * var_R**(6/5))**(5/6))) - 1
BER = np.zeros((len(L), len(Ps)))

for i in range(len(L)):
    for j in range(len(Ps)):
        for k in range(len(w20)):
            term = (m_mod / 4) * Ps[j] * a[i] * np.exp(np.sqrt(2 * var_s[i]) * x20[k] - var_s[i] / 2)
            term = term / np.max(term)  # normalize the term
            var_N_k = 2 * q_charge * (gain_apd**2) * FA * R * delta_f * term + (4 * Boltzmann * T_Rx * Fn * delta_f / R_load)
            BER[i, j] = BER[i, j] + w20[k] * 0.5 * erfc((m_mod * R * gain_apd * Ps[j] * a[i]) / (4 * np.sqrt(var_N_k)) * np.exp(np.sqrt(2 * var_s[i]) * x20[k] - var_s[i] / 2))
        BER[i, j] = BER[i, j] * (1 / np.sqrt(np.pi))

plt.figure(figsize=(10, 6))
for i in range(len(L)):
    plt.semilogy(Ps_dB, BER[i, :], label=f"L = {L[i]}m")

plt.xlabel("Transmitted Power (dB)")
plt.ylabel("BER")
plt.xlim([0, 50])
plt.ylim([1e-10, 1])  # set y-axis limits
plt.title("BER vs Transmitted Power for varying Distance")
plt.legend()
plt.grid(True)
plt.show()