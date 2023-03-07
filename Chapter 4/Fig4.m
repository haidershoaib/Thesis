clc;
clear all;
close all;

V_max = 30; % m/s
R_th =  5*10^7; % B/sec
W_c = 40*10^6; % Hz
P_R = 1; % Watts
mu_max = 0.01; % Max BS density
h_d = 3/(mu_max*V_max); % Normalized handover rate

% SINR parameters
% alpha = [2, 4]; % path loss exponent
alpha = 3;
alpha1 = 2;
G_tx = 1; % Gain at transmitter
G_rx = 1; % Gain at reciever
c = 3*10^8; % m/s
f_R = 2.1*10^9; % Hz
NP=W_c*273*1.38*(10)^-23; % watts/m^2
gamma_R = G_tx*G_rx*(c/(4*pi*f_R))^2;
lambda_i = 1;
lambda_s = 1;
Pj = 1;

% gamma_bar = (gamma_R*P_R)/NP;
mu = 0.001:0.0001:0.01; % BS densities 
L = 2000; % Length of road
d_safe = 5; % BS safety distance
h_bs = 8; % BS height
for i=1:length(mu)
    N = mu(i)*L;
    for j = 1:N-1
        d(j) = (h_bs^2 + d_safe^2 + ((2*j+1)^2/(4*mu(i)^2)))^(-alpha/2);
        d1(j) = (h_bs^2 + d_safe^2 + ((2*j+1)^2/(4*mu(i)^2)))^(-alpha1/2);
    end
    lambda = lambda_s./(gamma_R.*d);

    rho_0 = ((1./(2.*mu(i))).^2 + h_bs^2 + d_safe^2)^(-alpha/2);
    rho_1 = ((1./(2.*mu(i))).^2 + h_bs^2 + d_safe^2)^(-alpha1/2);
    lambda2 = lambda_s./(Pj*gamma_R.*rho_0);
    lambda21 = lambda_s./(Pj*gamma_R.*rho_1);

    fun = @(z)(exp(-NP*z).*(M_I(z, mu(i), L, Pj*gamma_R, lambda_i, d_safe, h_bs, alpha).*(1-M_S(z, mu(i), d_safe, h_bs, gamma_R, lambda_s, alpha)))./(z));
    R_avg(i) = integral(fun, 0, inf);

    fun1 = @(z)(exp(-NP*z).*(M_I(z, mu(i), L, Pj*gamma_R, lambda_i, d_safe, h_bs, alpha1).*(1-M_S(z, mu(i), d_safe, h_bs, gamma_R, lambda_s, alpha1)))./(z));
    R_avg1(i) = integral(fun1, 0, inf);


%     fun2 = @(z)((1-hypo(z, lambda, lambda21))./(1+z));
%     erg_cap2(i) = W_c.*(1/(log(2))).*integral(fun2, 0, inf);
end

for j=1:length(mu)
    if (W_c*R_avg(j) < R_th)
        break;
    else
      V_data(j) = (1./(h_d.*mu(j))).*(1-R_th./(W_c.*R_avg(j)));
      V_data1(j) = (1./(h_d.*mu(j))).*(1-R_th./(W_c.*R_avg1(j)));
    end
end

tau = 0.003;
epsilon = 0.01;
sigma_LN = 1;
mu_LN = 0;
V_safe = (1/tau)*exp(sigma_LN*sqrt(2)*erfinv(2*epsilon -1) + mu_LN);
vmin = min(V_max, V_safe);

tau1 = 0.004;
V_safe1 = (1/tau1)*exp(sigma_LN*sqrt(2)*erfinv(2*epsilon -1) + mu_LN);
vmin1 = min(V_max, V_safe1);

subplot(2,1,1)
plot(mu, V_data, 'k', 'LineWidth', 1.2)
hold on;
% plot(mu, V_data1, 'k', 'LineWidth', 1.2)
% hold on
plot(mu, vmin+mu*0, 'r', 'LineWidth', 1.2);
xlabel('BS Density (\mu) [BSs/m]')
ylabel('V_{data} [m/s]')
grid on;
ylim([0 35])
legend('V_{data} for \alpha = 3', 'V_{min} = min(V_{safe}, V_{max}) for \tau = 0.003 s');

subplot(2,1,2)
plot(mu, V_data1, 'k', 'LineWidth', 1.2)
hold on;
% plot(mu, V_data1, 'k', 'LineWidth', 1.2)
% hold on
plot(mu, vmin1+mu*0, 'r', 'LineWidth', 1.2);
xlabel('BS Density (\mu) [BSs/m]')
ylabel('V_{data} [m/s]')
grid on;
legend('V_{data} for \alpha = 4', 'V_{min} = min(V_{safe}, V_{max}) for \tau = 0.004 s');
% ylim([0 35])
% plot(mu, V_data1)

function h = hypo(z, lambda_i, lambda_s)
cdf = 0;
for i = 1:length(lambda_i)
    W = 1;
    for j = 1:length(lambda_i)
        if j == i
            continue;
        end
        W = W*(lambda_i(j)/(lambda_i(j)-lambda_i(i)));
    end
    cdf = cdf + W.*(lambda_i(i)./(lambda_s.*z+lambda_i(i)));
end
h = 1 - cdf;
end

function m = M_I(z, mu, L, gamma_R, lambda_i, x, zh, alpha)
    m = 1;
    N = mu*L;
    for i=1:N-1
        d = (x^2 + zh^2 + ((2*i+1)^2/(4*mu^2)))^(-alpha/2);
        m = m.*(1./(1+z.*(((gamma_R.*d)./(lambda_i)))));
    end
end

function ms = M_S(z, mu, x, zh, gamma_R, lambda_s, alpha)
    rho_0 = sqrt((1./(2.*mu)).^2 + zh^2 + x^2);
    ms = 1./(1+z.*(((gamma_R*rho_0.^(-alpha))./(lambda_s))));
end