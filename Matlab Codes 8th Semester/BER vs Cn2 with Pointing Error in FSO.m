% Constants
q_charge = 1.602e-19;  % C
m_mod = 1;
R_load = 1e3;  % Ohm
Fn = 2;
Rb = [1e9];
kA = 0.7;
R = 1;
D = 0.02;  % m
beta_v = (0.1 * log(10)) / 10000;  % dB/km ---> 1/km
phi_divergence = 1e-3;  % rad
wav_lambda = 1.55e-6;  % m
T_Rx = 300;  % K
Cn2 = linspace(0.5e-14, 4e-14, 100); % 20 values in between 0.1e-14 and 4e-14 on linear scale gap
theta = 1e-6;
sigma = 1e-6;
pointing_error_factor = exp(-(theta^2) / (2 * sigma^2));
Ps = 10^(-3);  % W
L = 1e3;  % m
w20 = [2.22939364554e-13, 4.39934099226e-10, 1.08606937077e-7, 7.8025564785e-6, 0.000228338636017, 0.00324377334224, 0.0248105208875, 0.10901720602, 0.286675505363, 0.462243669601, 0.462243669601, 0.286675505363, 0.10901720602, 0.0248105208875, 0.00324377334224, 0.000228338636017, 7.8025564785e-6, 1.08606937077e-7, 4.39934099226e-10, 2.22939364554e-13];
x20 = linspace(-5.38748089001, 5.38748089001, 20);
gain_apd = 15;
delta_f = Rb / 2;
a = ((D / (phi_divergence * L))^2) * exp(-beta_v * L);

FA = kA * gain_apd + (1 - kA) * (2 - (1 / gain_apd));
var_R = 1.23 * Cn2 * ((2 * pi / wav_lambda)^(7/6)) * (L^(11/6));
d = sqrt(((2 * pi / wav_lambda) * D^2) / (4 * L));
var_s = exp((0.49 * var_R) / ((1 + 0.18 * d^2 + 0.56 * var_R^(6/5))^(7/6)) + (0.51 * var_R) / ((1 + 0.9 * d^2 + 0.62 * (d^2) * var_R^(6/5))^(5/6))) - 1;
BER = zeros(length(Rb), length(Cn2));
Legend = cell(1,2);

for i = 1:length(Rb)
    for j = 1:length(Cn2)
        for k = 1:length(w20)
            var_N_k = 2 * q_charge * (gain_apd^2) * FA * R * delta_f(i) * (m_mod / 4) * Ps * a * exp(sqrt(2 * var_s(j)) * x20(k) - var_s(j) / 2) + (4 * 1.380649e-23 * T_Rx * Fn * delta_f(i) / R_load);  % Boltzmann constant in J/K
            BER(i, j) = BER(i, j) + w20(k) * 0.5 * erfc((m_mod * R * gain_apd * Ps * a) / (4 * sqrt(var_N_k)) * exp(sqrt(2 * var_s(j)) * x20(k) - var_s(j) / 2));
        end
        BER(i, j) = BER(i, j) * (1 / sqrt(pi));
    end
    semilogy(Cn2, BER(i, :), 'Color', 'y');
    hold on;
    Legend{1} = "Without Pointing Error";
    Ps = Ps*pointing_error_factor;
    for j = 1:length(Cn2)
        for k = 1:length(w20)
            var_N_k = 2 * q_charge * (gain_apd^2) * FA * R * delta_f(i) * (m_mod / 4) * Ps * a * exp(sqrt(2 * var_s(j)) * x20(k) - var_s(j) / 2) + (4 * 1.380649e-23 * T_Rx * Fn * delta_f(i) / R_load);  % Boltzmann constant in J/K
            BER(i, j) = BER(i, j) + w20(k) * 0.5 * erfc((m_mod * R * gain_apd * Ps * a) / (4 * sqrt(var_N_k)) * exp(sqrt(2 * var_s(j)) * x20(k) - var_s(j) / 2));
        end
        BER(i, j) = BER(i, j) * (1 / sqrt(pi));
    end
    semilogy(Cn2, BER(i, :), 'Color', 'b');
    Legend{2} = "With Pointing Error";
end

xlabel('Cn2');
ylabel('BER');
legend(Legend);
grid on;