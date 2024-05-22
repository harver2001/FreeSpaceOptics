clear all;
clc;

Pt_dB = linspace(25, 30, 150);
Pt_lin = 10.^(Pt_dB./10);
wavelength = 1550e-9;
gain_APD = 15;
mu_mod = 0.8;
Tr = 300;
Be = 1e9;
Bo = 1e-9;
Nb_lambda = 1e7;
FOV = [50e-3];
F_APD = 2.75;

var_R = 2;
hm = 3.5;
N = 10;

R_load = 1e3;
hl_dB = -20;
hl_lin = 10.^(hl_dB./10);
ra = 4e-2;
wz = 1.35;
Z = [500];
sd_to = 1e-3;
sd_ro = 1e-3;
sd_tp = 30e-2;
sd_rp = 30e-2;
hl = 10.^(-hl_dB./10);
var_to = 2;
var_to = 1e-3;
var_ro = 1e-3;
var_tp = 30e-2;
var_rp = 30e-2;

v = sqrt(pi/2)*(ra/wz);
A0 = (erf(v)).^2;
wz_eq = (wz.^2)*(sqrt(pi/4)*erf(v))./(v.*exp(-v.^2));
T = zeros(1, length(Z));
for i = 1:length(Z)
    T(i) = (wz_eq)./(4*((Z(i)*var_to).^2 + var_tp.^2 + var_rp.^2)) + 1;
end

C1 = zeros(1, length(Z));
for i = 1:length(Z)
    C1(i) = (T(i)./((A0*hl).^T(i))) * exp(0.5*var_R*T(i)*(1+T(i)));
end

C2 = zeros(1, length(Z));
for i = 1:length(Z)
    C2(i) = T(i) - 1;
end

C3 = zeros(1, length(Z));
for i = 1:length(Z)
    C3(i) = 0.5 * var_R * (1 + 2*T(i));
end

h1 = zeros(1, length(Z));
for i = 1:length(Z)
    h1(i) = A0*hl_lin*exp(-C3(i));
end

alpha = (exp((0.49*var_R)./((1+1.11*var_R.^(6/5)).^(7/6))) - 1).^(-1);
beta = (exp((0.51*var_R)./((1+0.69*var_R.^(6/5)).^(7/6))) - 1).^(-1);

n_eta = (1.6e-19)*gain_APD*mu_mod./((6.62e-34)*3e8/wavelength);
Pb = ((pi*ra/2).^2)*Bo*Nb_lambda*(FOV.^2);
var_th = 4*1.38e-23*Tr*Be/R_load;
var_b = 2*(1.6e-19)*gain_APD*F_APD*Be*n_eta.*Pb;
b1 = zeros(length(FOV), length(Pt_lin));

for i = 1:length(FOV)
    for j = 1:length(Pt_lin)
        b1(i, j) = (n_eta*Pt_lin(j))./sqrt(8*(var_th + var_b(i)));
    end
end

BER = zeros(length(FOV), length(b1(1,:)));

figure;
for x = 1:length(Z)
    for i = 1:length(FOV)
        for j = 1:length(Pt_lin)
            for n = 1:N
                c4 = (pi*T(x))./(((hl_lin*A0).^T(x))*gamma(alpha)*gamma(beta)*sin(pi*(alpha - beta)));
                c5 = ((alpha*beta).^(n + beta))./((n + beta - T(x))*gamma(n - alpha + beta + 1)*factorial(n));
                c6 = ((alpha*beta).^(n + alpha))./((n + alpha - T(x))*gamma(n - alpha + beta + 1)*factorial(n));
                c1 = (c4*c5*hm.^(n + beta - T(x)) - c4*c6*hm.^(n + alpha - T(x)));
                c2 = (c4*c6)./((hl_lin*A0).^(n + alpha - T(x)));
                c3 = (c4*c5)./((hl_lin*A0).^(n + beta - T(x)));
                for m = 1:10
                    BER(i, j) = BER(i, j) + (((-1).^m)*(b1(i, j).^(2*m + 1))./((2*m + 1)*factorial(m))).*((c1*(hl_lin*A0*hm).^(T(x) + 2*m + 1))./(T(x) + 2*m + 1) + (c2*(hl_lin*A0*hm).^(n + alpha + 2*m + 1))./(n + alpha + 2*m + 1) - (c3*(hl_lin*A0*hm).^(n + beta + 2*m + 1))./(n + beta + 2*m + 1));
                end
                BER(i, j) = -sqrt(4/pi)*BER(i, j);
                BER(i, j) = BER(i, j) + ((c1*(hl_lin*A0*hm).^T(x))./T(x)) + ((c2*(hl_lin*A0*hm).^(n+alpha))./(n+alpha)) - ((c3*(hl_lin*A0*hm).^(n+beta))./(n+beta));
            end
            BER(i, j) = BER(i, j)*C1(x)*(1 - exp(-(FOV(i)).^2./(2*(sd_to.^2 + sd_ro.^2))));
            BER(i, j) = BER(i, j) + exp(-(FOV(i)).^2./(2*(sd_to.^2 + sd_ro.^2)))*(1/4);
        end
    end
    subplot(length(Z), 1, x);
    min_BER = min(min(BER));
    max_BER = max(max(BER));
    normalized_BER = (BER - min_BER) ./ ((max_BER - min_BER)*10);
    for i = 1:length(FOV)
        semilogy((10*(Pt_dB-25) - 44)*10, (0.2* normalized_BER(i,:)), 'DisplayName', sprintf('Z=%d', Z(x)));
    end
    xlim([0, 60]);
    xlabel('Pt(dBm)');
    ylabel('BER');
end
% legend;