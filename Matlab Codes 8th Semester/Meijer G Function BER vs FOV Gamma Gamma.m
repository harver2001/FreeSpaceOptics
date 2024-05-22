alpha = 1.3;
beta = 3.2;
g = 1;
omega_prime = 1;
sigma_a = 1;
Am = 1;
Variable1 = [1, 1];

FOV = linspace(0.005,0.05,100);

Pe_const = [];
k = 1;
for j = 1:length(FOV)
    for i = 1:round(beta)-1
        am = (nchoosek(round(beta)-1,i)*(g*beta+omega_prime)^(1-((i+1)/2))*(omega_prime/g)^(i) *(alpha/beta)^((i+1)/2))/(factorial(i));
        Pe_const_temp = (am*(alpha*beta/(g*beta + omega_prime))^(-(alpha+i+1)/2))*meijerG([[], [1, 0.5]], [[0], Variable1], 1.0);
    end
    Pe_const = [Pe_const, (0.5*(exp(-(FOV(j)^2)/(2*sigma_a^2))) + (1-exp(-(FOV(j)^2)/(2*sigma_a^2)))*(Variable1(1)*Am/(4*sqrt(pi)))) / k];
    k = k + 7;
end

semilogy(FOV * 1000, Pe_const);
xlabel('Theta FOV (mrad)');
ylabel('BER');
title('Meijer Function');