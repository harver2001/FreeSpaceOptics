alpha = 1.3;
beta = 3.2;
g = 1;
omega_prime = 1;
sigma_a = 1;
Am = 1;
Variable1 = [1, 1];

SNR_dB = linspace(0, 30, 100);  
SNR = 10 .^ (SNR_dB / 10);  

Pe_const = [];
for snr = SNR
    Pe_const_temp = 0;
    for i = 1:round(beta)-1
        am = (nchoosek(round(beta)-1, i) * (g*beta + omega_prime) ^ (1 - ((i+1)/2)) * (omega_prime/g) ^ i * (alpha/beta) ^ ((i+1)/2)) / factorial(i);
        Pe_const_temp = Pe_const_temp + (am * (alpha*beta / (g*beta + omega_prime)) ^ (-(alpha+i+1)/2)) * meijerG([[], [1, 0.5]], [[0], Variable1], 1.0);
    end
    Pe_const = [Pe_const, (0.5 * exp(-snr/10) + (1 - exp(-snr/10)) * (Variable1(1) * Am / (4 * sqrt(pi)))) / 10^3];
end

semilogy(SNR_dB, Pe_const);
xlabel('SNR (dB)');
ylabel('BER');
xlim([0,16]);
title('BER vs. SNR');