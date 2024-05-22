clear all;
close all;

q_charge = 1.602e-19;  % C
gain_apd_values = [15];
f_apd = 2.75;

h = linspace(0,1e-3,50);
ra = 4e-2;
wz = 1;
var_R = 0.2;
var_to = 1e-3;
var_ro = 1e-3;
var_tp = 30e-2;
var_rp = 30e-2;
hl_dB = 10;
hl = 10^(-hl_dB/10);
Z = 500;
Z1 = [500, 1000, 1500, 2000, 2500, 3000];
FOV = linspace(0.005,0.05,100);
Theta_FOV = 5e-3;
v = sqrt(pi/2)*(ra/wz);
wz_eq = wz^2*(sqrt(pi)*erf(v))/(2*v*exp(-v^2));
A0 = (erf(v))^2;
T = (wz_eq)/(4*((Z*var_to)^2 + var_tp^2 + var_rp^2)) + 1;
C1 = (T/((A0*hl)*T)) * exp(0.5*var_R*T*(1+T));
C2 = T - 1;
C3 = 0.5 * var_R * (1 + 2*T);
a=[5/24,4/24,1/24];
az = [2,11/20,1/2];
Rl = 1e3; % Load Resistance
Tr = 300; % Receiver Temperature
lamda = 1550e-9;
c = 3e8;
v1 = c/lamda;
G = 30; % Average APD Gain
Mew = 0.8;
Be = 1e9; % electrical bandwidth
B0 = 1e-9; % optical filter bandwidth
Nb = 1e-3; % spectral radiance of background radiation 
eta = q_charge*G*Mew/(Planck*v1);

Pb = (pi^2 * ra^2 * B0 * Nb * FOV.^2) / 4;
var_th = (4*Boltzmann*Tr*(1e9)) / (Rl);
var_bg = 2*q_charge*G*f_apd*Be*eta*Pb;
var_n = var_bg + var_th;
C4 = [];

for i = 1:length(az)
    x = -az(i)*log(A0*hl)+az(i)*C3 - var_R*(C2+1)/2;
    C4 = [C4, x];
end

Pt = 1:5:30;
YTh = 1e6;
R = 1;
Hth = zeros(1, length(Pt));

for i = 1:length(Pt)
    Hth(i) = (i+1.5)*1e-7;
end

figure;
for h = 1:length(Hth)
    Pout = [];
    for i = 1:length(FOV)
        temp = 0;
        for j = 1:length(a)
            temp = temp + (a(j)/2)*sqrt(pi/az(j))*exp((C4(j)^2 - (az(j)*(log(A0*hl) - C3))^2)/az(j))*(1 + erf(sqrt(az(j))*log(Hth(h)) + C4(j)/sqrt(az(j))));
        end
        temp = (Hth(h)^(C2 + 1))/(C2 + 1) - temp;
        temp = temp*C1*(1 - exp(-(FOV(i)^2)/(2*(var_to^2 + var_ro^2)))) + exp(-(FOV(i)^2)/(2*(var_to^2 + var_ro^2)));
        Pout = [Pout, temp];
    end
    semilogy(FOV*1000, Pout, 'DisplayName', ['Pt=', num2str(Pt(length(Pt)-h+1)), ' dBm']);
    hold on;
end
xlabel('FOV(mrad)');
ylabel('Outage Probablity');
xlim([7,12]);
ylim([10^-7.5,10^-5]);
legend();
hold off;