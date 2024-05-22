import numpy as np
import matplotlib.pyplot as plt
from scipy.special import comb
from mpmath import meijerg
import math

alpha = 1.3
beta = 3.2
g = 1
omega_prime = 1
sigma_a = 1
Am = 1
Variable1 = [1, 1]

FOV = np.linspace(0.005,0.05,100)

Pe_const = []
k = 1
for j in range(len(FOV)):
    for i in range(1,round(beta)):
        am = (comb(round(beta)-1,i)*(g*beta+omega_prime)**(1-((i+1)/2))*(omega_prime/g)**(i) *(alpha/beta)**((i+1)/2))/(math.factorial(i))
        Pe_const_temp = (am*(alpha*beta/(g*beta + omega_prime))**(-(alpha+i+1)/2))*meijerg([[], [1, 0.5]], [[0], Variable1], 1.0)
    Pe_const.append((0.5*(np.exp(-(FOV[j]**2)/(2*sigma_a**2))) + (1-np.exp(-(FOV[j]**2)/(2*sigma_a**2)))*(Variable1[0]*Am/(4*np.sqrt(np.pi)))) / k)
    k += 7

plt.semilogy(FOV * 1000, Pe_const)
plt.xlabel('Theta FOV (mrad)')
plt.ylabel('BER')
plt.title('Meijer Function')
plt.show()