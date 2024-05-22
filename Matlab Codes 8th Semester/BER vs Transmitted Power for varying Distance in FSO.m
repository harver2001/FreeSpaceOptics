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
Cn2 = 6e-14;
Ps = logspace(-1, 7, 20);  % W
Ps_dB = 10 * log10(Ps);
L = [1e3, 1.3e3, 1.6e3];  % m
w20 = [2.22939364554e-13, 4.39934099226e-10, 1.08606937077e-7, 7.8025564785e-6, 0.000228338636017, 0.00324377334224, 0.0248105208875, 0.10901720602, 0.286675505363, 0.462243669601, 0.462243669601, 0.286675505363, 0.10901720602, 0.0248105208875, 0.00324377334224, 0.000228338636017, 7.8025564785e-6, 1.08606937077e-7, 4.39934099226e-10, 2.22939364554e-13];
x20 = linspace(-5.38748089001, 5.38748089001, 20);
gain_apd = 15;
delta_f = Rb / 2;
a = ((D ./ (phi_divergence * L)).^2) .* exp(-beta_v * L);

FA = kA * gain_apd + (1 - kA) * (2 - (1 / gain_apd));
var_R = 1.23 * Cn2 * ((2 * pi / wav_lambda)^(7/6)) * (L.^(11/6));
d = sqrt(((2 * pi / wav_lambda) * D^2) / (4 * L));
var_s = exp((0.49 * var_R) ./ ((1 + 0.18 * d.^2 + 0.56 * var_R.^(6/5)).^(7/6)) + (0.51 * var_R) ./ ((1 + 0.9 * d.^2 + 0.62 * (d.^2) * var_R.^(6/5)).^(5/6))) - 1;
BER = zeros(length(L), length(Ps));

for i = 1:length(L)
    for j = 1:length(Ps)
        for k = 1:length(w20)
            term = (m_mod / 4) * Ps(j) * a(i) * exp(sqrt(2 * var_s(i)) * x20(k) - var_s(i) / 2);
            term = term / max(term);  % normalize the term
            var_N_k = 2 * q_charge * (gain_apd^2) * FA * R * delta_f * term + (4 * 1.38064852e-23 * T_Rx * Fn * delta_f / R_load);
            BER(i, j) = BER(i, j) + w20(k) * 0.5 * erfc((m_mod * R * gain_apd * Ps(j) * a(i)) / (4 * sqrt(var_N_k)) * exp(sqrt(2 * var_s(i)) * x20(k) - var_s(i) / 2));
        end
        BER(i, j) = BER(i, j) * (1 / sqrt(pi));
    end
end

figure('Position', [10 10 900 600])
for i = 1:length(L)
    semilogy(Ps_dB, BER(i, :), 'DisplayName', sprintf('L = %dm', L(i)));
    hold on;
end

xlabel('Transmitted Power (dB)')
ylabel('BER')
xlim([0, 50])
ylim([1e-10, 1])
title('BER vs Transmitted Power for varying Distance')
legend()
grid on
hold off