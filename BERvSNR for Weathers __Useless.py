import numpy as np
from scipy.constants import Boltzmann
from scipy.special import erfc
import matplotlib.pyplot as plt

# I_sun and I_sky
N_lambda = 1e-3*(1e4*1e6) #W/cm2µmSr --> W/m3Sr
W_lambda = 0.055*(1e4*1e6) #W/cm2µm --> W/m3
del_lambda_OBF = 1e-9 #nm --> m
FOV = 0.6 #rad

Isky = N_lambda*del_lambda_OBF*np.pi*(FOV**2)/4
Isun = W_lambda*del_lambda_OBF

# var_bg and var_th
q_charge = 1.602e-19 #C
B_ef = 155e6 #Hz
R = 1
Te = 300 #Kelvin
R_load = 50 #Ohm

var_bg = 2*q_charge*B_ef*R*(Isky+Isun)
var_th = (4*Boltzmann*Te*B_ef)/R_load

# K0 and var_ryotov
mod_index = 1
A = 1
Cn2 = 0.75e-14 #m(-2/3)
Io = np.logspace(-10,-4,30)
Pd = (A**2)/2
lambda_ = 850e-9
L = 1e3 #m

K0 = np.zeros(Io.shape)
var_ryotov = 1.23*Cn2*((2*np.pi/lambda_)**(7/6))*(L**(11/6))

# Fog and Visibility effects
lambda0 = 550e-9 #m
Visibility = np.linspace(1,100,15)*1e3 #km --> m
attenuation_coeff = np.where(Visibility < 6, (10*np.log10(5)*(lambda_/lambda0)**(-0.585*Visibility**(1/3)))/Visibility, np.where(Visibility < 50, (10*np.log10(5)*(lambda_/lambda0)**(-1.3))/Visibility, (10*np.log10(5)*(lambda_/lambda0)**(-1.6))/Visibility))

# w20 and x20 vectors
w20 = np.array([2.22939364554e-13,4.39934099226e-10,1.08606937077e-7, 7.8025564785e-6,0.000228338636017,0.00324377334224,0.0248105208875, 0.10901720602,0.286675505363,0.462243669601, 0.462243669601,0.286675505363,0.10901720602,0.0248105208875, 0.00324377334224,0.000228338636017,7.8025564785e-6,1.08606937077e-7, 4.39934099226e-10,2.22939364554e-13])
x20 = np.array([-5.38748089001,-4.60368244955,-3.94476404012,-3.34785456738,-2.78880605843,-2.25497400209,-1.73853771212,-1.2340762154,-0.737473728545,-0.245340708301,0.245340708301,0.737473728545,1.2340762154,1.73853771212,2.25497400209,2.78880605843,3.34785456738,3.94476404012,4.60368244955,5.38748089001])

# BER calculation
BER = np.zeros((len(attenuation_coeff),len(Io)))
SNR = R_load*((R*Io)**2)/(4*Boltzmann*Te*B_ef)
SNRdB = 10*np.log10(SNR)
Legend = [None]*len(attenuation_coeff)

plt.figure()
for k in range(len(attenuation_coeff)):
    SNR = R_load*((R*Io*np.exp(-attenuation_coeff[k]*L))**2)/(4*Boltzmann*Te*B_ef)
    SNRdB = 10*np.log10(SNR)
    K0 = (((mod_index*R*Io*np.exp(-attenuation_coeff[k]*L))**2)*Pd)/(var_bg + var_th)
    for i in range(len(Io)):
        for j in range(len(x20)):
            BER[k,i] = BER[k,i] + w20[j]*0.5*erfc(np.sqrt(K0[i])*np.exp(x20[j]*np.sqrt(2)*np.sqrt(var_ryotov) - var_ryotov/2)/np.sqrt(2))
    plt.semilogy(SNRdB, BER[k,:])
    # Legend[k] = "Vis. = " + str(Visibility[k])

plt.show()