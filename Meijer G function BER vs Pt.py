import numpy as np
import matplotlib.pyplot as plt
from scipy.special import comb
from mpmath import meijerg
import math

# Constants and Parameters
alpha = 1.3
beta = 3.2
g = 1
omega_prime = 1
sigma_a = 1
Am = 1
Variable1 = [1, 1]

# Simulate SNR values
SNR_dB = np.linspace(0, 30, 100)  # SNR range from 0 dB to 30 dB
SNR = 10 ** (SNR_dB / 10)  # Convert SNR from dB to linear scale

Pe_const = []
for snr in SNR:
    Pe_const_temp = 0
    for i in range(1, round(beta)):
        am = (comb(round(beta)-1, i) * (g*beta + omega_prime) ** (1 - ((i+1)/2)) * (omega_prime/g) ** i * (alpha/beta) ** ((i+1)/2)) / math.factorial(i)
        Pe_const_temp += (am * (alpha*beta / (g*beta + omega_prime)) ** (-(alpha+i+1)/2)) * meijerg([[], [1, 0.5]], [[0], Variable1], 1.0)
    Pe_const.append((0.5 * np.exp(-snr/10) + (1 - np.exp(-snr/10)) * (Variable1[0] * Am / (4 * np.sqrt(np.pi)))) / 10**3)

plt.semilogy(SNR_dB, Pe_const)
plt.xlabel('SNR (dB)')
plt.ylabel('BER')
plt.xlim([0,16])
plt.title('BER vs. SNR')
plt.show()