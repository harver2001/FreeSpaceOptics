import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc

# Constants
q_charge = 1.602e-19  # C
m_mod = 1
R_load = 1e3  # Ohm
Fn = 2
# Rb = np.array([2e9, 1e9, 5e8, 2e8, 1e8])
Rb = np.array([1e9])
kA = 0.7
R = 1
D = 0.02  # m
beta_v = (0.1 * np.log(10)) / 10000  # dB/km ---> 1/km
phi_divergence = 1e-3  # rad
wav_lambda = 1.55e-6  # m
T_Rx = 300  # K
Cn2 = np.linspace(0.5e-14, 4e-14, 100) # 20 values in between 0.1e-14 and 4e-14 on linear scale gap
theta = 1e-6
sigma = 1e-6
pointing_error_factor = np.exp(-(theta**2) / (2 * sigma**2))
Ps = 10**(-3)  # W
# Ps = 1e-3
L = 1e3  # m
w20 = np.array([2.22939364554e-13, 4.39934099226e-10, 1.08606937077e-7, 7.8025564785e-6, 0.000228338636017, 0.00324377334224, 0.0248105208875, 0.10901720602, 0.286675505363, 0.462243669601, 0.462243669601, 0.286675505363, 0.10901720602, 0.0248105208875, 0.00324377334224, 0.000228338636017, 7.8025564785e-6, 1.08606937077e-7, 4.39934099226e-10, 2.22939364554e-13])
x20 = np.linspace(-5.38748089001, 5.38748089001, 20)
gain_apd = 15
delta_f = Rb / 2
a = ((D / (phi_divergence * L))**2) * np.exp(-beta_v * L)

FA = kA * gain_apd + (1 - kA) * (2 - (1 / gain_apd))
var_R = 1.23 * Cn2 * ((2 * np.pi / wav_lambda)**(7/6)) * (L**(11/6))
d = np.sqrt(((2 * np.pi / wav_lambda) * D**2) / (4 * L))
var_s = np.exp((0.49 * var_R) / ((1 + 0.18 * d**2 + 0.56 * var_R**(6/5))**(7/6)) + (0.51 * var_R) / ((1 + 0.9 * d**2 + 0.62 * (d**2) * var_R**(6/5))**(5/6))) - 1
BER = np.zeros((len(Rb), len(Cn2)))
Legend = [None]*2

for i in range(len(Rb)):
    for j in range(len(Cn2)):
        for k in range(len(w20)):
            var_N_k = 2 * q_charge * (gain_apd**2) * FA * R * delta_f[i] * (m_mod / 4) * Ps * a * np.exp(np.sqrt(2 * var_s[j]) * x20[k] - var_s[j] / 2) + (4 * 1.380649e-23 * T_Rx * Fn * delta_f[i] / R_load)  # Boltzmann constant in J/K
            BER[i, j] = BER[i, j] + w20[k] * 0.5 * erfc((m_mod * R * gain_apd * Ps * a) / (4 * np.sqrt(var_N_k)) * np.exp(np.sqrt(2 * var_s[j]) * x20[k] - var_s[j] / 2))
        BER[i, j] = BER[i, j] * (1 / np.sqrt(np.pi))
    # plt.subplot(len(Rb), 1, i + 1)
    plt.semilogy(Cn2, BER[i, :], color = 'y')
    Legend[0] = "Without Pointing Error"
    Ps = Ps*pointing_error_factor
    for j in range(len(Cn2)):
        for k in range(len(w20)):
            var_N_k = 2 * q_charge * (gain_apd**2) * FA * R * delta_f[i] * (m_mod / 4) * Ps * a * np.exp(np.sqrt(2 * var_s[j]) * x20[k] - var_s[j] / 2) + (4 * 1.380649e-23 * T_Rx * Fn * delta_f[i] / R_load)  # Boltzmann constant in J/K
            BER[i, j] = BER[i, j] + w20[k] * 0.5 * erfc((m_mod * R * gain_apd * Ps * a) / (4 * np.sqrt(var_N_k)) * np.exp(np.sqrt(2 * var_s[j]) * x20[k] - var_s[j] / 2))
        BER[i, j] = BER[i, j] * (1 / np.sqrt(np.pi))
    # plt.title(f"Rb = {Rb[i]}")
    plt.semilogy(Cn2, BER[i, :], color = 'b')
    Legend[1] = "With Pointing Error"

plt.xlabel('Cn2')
plt.ylabel('BER')
plt.legend(Legend)
plt.show()