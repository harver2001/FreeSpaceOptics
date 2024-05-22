import numpy as np
import matplotlib.pyplot as plt
import math

var_th = 1
var_bg = 1
Yth = 1e6
R = 1
C4 = []
az = np.array([1, 2, 3])  # Assuming az is an array of length 3
A0 = 1  # Assuming A0 is 1
hl = 1  # Assuming hl is 1
var_R = 1  # Assuming var_R is 1
C2 = 1  # Assuming C2 is 1
Pt = np.array([0.005,0.010,0.015,0.020,0.025,0.030,0.035,0.040])
eta = 1  # Assuming eta is 1
var_th = 1  # Assuming var_th is 1
var_bg = 1  # Assuming var_bg is 1
FOV = np.array([1, 2, 3])  # Assuming FOV is an array of length 3
var_to = 1  # Assuming var_to is 1
var_ro = 1  # Assuming var_ro is 1
wz_eq = 1  # Assuming wz_eq is 1
var_tp = 1  # Assuming var_tp is 1
var_rp = 1  # Assuming var_rp is 1
Z = np.array([100, 200, 300, 400, 500])
BER = []

for i in range(len(az)):
    T = (wz_eq)/(4*((Z[i]*var_to)**2 + var_tp**2 + var_rp**2)) + 1
    C1 = (T/((A0*hl)*T)) * np.exp(0.5*var_R*T*(T+1))
    C2 = T - 1
    C3 = 0.5 * var_R * (1 + 2*T)
    x = -az[i]*np.log(A0*hl)+az[i]*C3 - var_R*(C2+1)/2
    C4.append(x)

b1 = (eta * Pt[:, np.newaxis]) / np.sqrt(8 * (var_th + var_bg))

plt.figure()  

for z in Z:
    T = (wz_eq)/(4*((z*var_to)**2 + var_tp**2 + var_rp**2)) + 1
    C1 = (T/((A0*hl)*T)) * np.exp(0.5*var_R*T*(T+1))
    C2 = T - 1
    C3 = 0.5 * var_R * (1 + 2*T)
    h1 = A0*hl*np.exp(-1*C3)
    Pe = [0] * len(Pt)  
    for i in range(len(Pt)):
        Pe[i] = 0.5 * np.exp((-1 * FOV[0] ** 2) / (2 * (var_to ** 2 + var_ro ** 2)))
        temp = 0
        for j in range(3):
            temp += 1  # Add your calculation here
        Pe[i] += temp
    BER.append(np.mean(Pe))

plt.plot(Z, BER) 
plt.xlabel('Z (Distance)')
plt.ylabel('BER')
# plt.ylim([10 ** -4, 0.11])
plt.show()