import numpy as np
from scipy.special import erfc

import matplotlib.pyplot as plt

W_lambda = 550e-4
N_lambda = 1e-3
delta_lambdaOBF = 1e-3
FOV_detector = 0.6

I_sun = W_lambda*delta_lambdaOBF
I_sky = N_lambda*delta_lambdaOBF*(4*FOV_detector**2)/np.pi

q = 1.602e-19
Responsivity = 1
B_ef = 155e6
k_boltzmann = 1.38e-23
Temperature = 500
B_ef = 155e6
R_load = 50
A = 1

var_BG = 2*q*Responsivity*(I_sky + I_sun)
var_Thermal = (4*k_boltzmann*Temperature*B_ef)/R_load
Pm = (A**2)/2

Irradiance = np.logspace(-10,-4,30)
IrradiancedBm = 10*np.log10(Irradiance)

modulation_index = 1

k1 = (((modulation_index*Responsivity*Irradiance)**2)*Pm)/(var_BG + var_Thermal)

Cn = 0.75e-14
lamda = 850e-9
L_link = np.linspace(100,5000,100)

var_ryotov = 1.23*Cn*((2*np.pi/lamda)**(7/6))*L_link**(11/6)

w20 = np.array([2.22939364554e-13,4.39934099226e-10,1.08606937077e-7, 7.8025564785e-6,0.000228338636017,0.00324377334224,0.0248105208875,0.10901720602,0.286675505363,0.462243669601, 0.462243669601,0.286675505363,0.10901720602,0.0248105208875,0.00324377334224,0.000228338636017,7.8025564785e-6,1.08606937077e-7,4.39934099226e-10,2.22939364554e-13])
x20 = np.array([-5.38748089001,-4.60368244955,-3.94476404012,-3.34785456738,-2.78880605843,-2.25497400209,-1.73853771212,-1.2340762154,-0.737473728545,-0.245340708301,0.245340708301,0.737473728545,1.2340762154,1.73853771212,2.25497400209,2.78880605843,3.34785456738,3.94476404012,4.60368244955,5.38748089001])
Legend = [None]*10

BER = np.zeros((100,30))
SNR = R_load*((Responsivity*Irradiance)**2)/(4*k_boltzmann*Temperature*B_ef)
SNRdB = 100*np.log10(SNR)

plt.figure()
selected_irradiance_indices = [20,21,22]  # Indices of the irradiance values we want to plot

for index in selected_irradiance_indices:
    j = index
    BER_j = []
    for i in range(len(L_link)):
        for k in range(len(w20)):
            BER[i,j] = BER[i,j] + w20[k]*0.5*erfc(np.sqrt(k1[j])*np.exp(x20[k]*np.sqrt(2)*np.sqrt(var_ryotov[i]) - var_ryotov[i]/2)/np.sqrt(2))
        BER_j.append(BER[i,j])
    plt.semilogy(L_link, BER_j, label=f"Irradiance {IrradiancedBm[j]} dBm")  # Plotting BER vs Link Distance

plt.xlabel('Link Distance')
plt.ylabel('BER')
plt.ylim([1e-10,0])
plt.xlim([0, 3000])
# plt.legend()
plt.show()