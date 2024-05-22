% Constants
q_charge = 1.602e-19;  % C
m_mod = 1;
R_load = 1e3;  % Ohm
Fn = 2;
Rb = 2e9;
kA = 0.7;
R = 1;
D = 0.02;  % m
beta_v = (0.1 * log(10)) / 10000;  % dB/km ---> 1/km
phi_divergence = 1e-3;  % rad
wav_lambda = 1.55e-6;  % m
T_Rx = 300;  % K
Cn2 = 6e-15;
Ps_values = logspace(-10, -0.1, 80); 
L = [1.6e3];  % m
x20 = linspace(-5.38748089001, 5.38748089001, 20);
w20 = [2.22939364554e-13, 4.39934099226e-10, 1.08606937077e-7, 7.8025564785e-6, 0.000228338636017, 0.00324377334224, 0.0248105208875, 0.10901720602, 0.286675505363, 0.462243669601, 0.462243669601, 0.286675505363, 0.10901720602, 0.0248105208875, 0.00324377334224, 0.000228338636017, 7.8025564785e-6, 1.08606937077e-7, 4.39934099226e-10, 2.22939364554e-13];
gain_apd_values = [15];
delta_f = Rb / 2;
a = ((D ./ (phi_divergence .* L)).^2) .* exp(-beta_v .* L); % attenuation factor

for gain_apd = gain_apd_values
    FA = kA * gain_apd + (1 - kA) * (2 - (1 / gain_apd));
    var_R = 1.23 * Cn2 * ((2 * pi / wav_lambda)^(7/6)) * (L^(11/6));
    d = sqrt(((2 * pi / wav_lambda) * D^2) / (4 * L));
    var_s = exp(var_R ^ 2);
    BER = zeros(length(L), length(Ps_values));
    SNR = zeros(length(L), length(Ps_values));

    for i = 1:length(L)
        for j = 1:length(Ps_values)
            Ps = Ps_values(j);
            SNR_PIN = R_load*((R*Ps*a(i))^2)/(4*Boltzmann*T_Rx*delta_f);  % Corrected line, moved inside the loop
            for k = 1:length(w20)
                var_N_k = 2 * q_charge * (gain_apd^2) * FA * R * delta_f * (m_mod / 4) * Ps * a(i) * exp(sqrt(2 * var_s(i)) * x20(k) - var_s(i) / 2) + (4 * Boltzmann * T_Rx * Fn * delta_f / R_load);
                BER(i, j) = BER(i, j) + w20(k) * 0.5 * erfc((m_mod * R * gain_apd * Ps * a(i)) / (4 * sqrt(var_N_k)) * exp(sqrt(2 * var_s(i)) * x20(k) - var_s(i) / 2));
                SNR(i, j) = (gain_apd * R_load * ((R * Ps * a(i) )^2)) / (2 * q_charge * (gain_apd * Fn * R * Ps * a(i) + 0) * delta_f + 4 * Boltzmann * T_Rx * delta_f / R_load);
            end
            BER(i, j) = BER(i, j) * (1 / sqrt(pi));
        end
        subplot(length(L), 1, i);
        semilogy(10 * log10(SNR(i, :)), BER(i, :)); hold on;
        legend('Without Pointing Error');
    end

    theta = 1e-6;
    sigma = 1e-6;
    pointing_error_factor = exp(-(theta^2) / (2 * sigma^2));
    Ps_values_adjusted = Ps_values * pointing_error_factor;  % Adjust Ps_values for pointing error
    for i = 1:length(L)
        for j = 1:length(Ps_values_adjusted)  % Use adjusted Ps_values here
            Ps = Ps_values_adjusted(j);
            SNR_PIN = R_load*((R*Ps*a(i))^2)/(4*Boltzmann*T_Rx*delta_f);  % Corrected line, moved inside the loop
            for k = 1:length(w20)
                var_N_k = 2 * q_charge * (gain_apd^2) * FA * R * delta_f * (m_mod / 4) * Ps * a(i) * exp(sqrt(2 * var_s(i)) * x20(k) - var_s(i) / 2) + (4 * Boltzmann * T_Rx * Fn * delta_f / R_load);
                BER(i, j) = BER(i, j) + w20(k) * 0.5 * erfc((m_mod * R * gain_apd * Ps * a(i)) / (4 * sqrt(var_N_k)) * exp(sqrt(2 * var_s(i)) * x20(k) - var_s(i) / 2));
                SNR(i, j) = (gain_apd * R_load * ((R * Ps * a(i) )^2)) / (2 * q_charge * (gain_apd * Fn * R * Ps * a(i) + 0) * delta_f + 4 * Boltzmann * T_Rx * delta_f / R_load);
            end
            BER(i, j) = BER(i, j) * (1 / sqrt(pi));
        end
        subplot(length(L), 1, i);
        semilogy(10 * log10(SNR(i, :)), BER(i, :)); hold on;
        legend('With Pointing Error');  
    end
end

xlabel('SNR (dB)');
ylabel('BER');
xlim([0,72]);
ylim([10^-5,1]);
grid on;
legend();  
hold off;