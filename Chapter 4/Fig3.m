clc;
clear all;
close all;

V_max = 30; % m/s
R_th =  5.9*10^7; % B/sec
W_c = 40*10^6; % Hz
P_R = 1; % Watts
mu_max = 0.02; % Max BS density
V = [0, 5, 15, 25];
h_d = 1/(mu_max*V_max); % Normalized handover rate
N_R = W_c*273*1.38*(10)^-23; % watts/m^2

% SINR parameters
% alpha = [2, 4]; % path loss exponent
alpha = 3;
G_tx = 1; % Gain at transmitter
G_rx = 1; % Gain at reciever
c = 3*10^8; % m/s
f_R = 2.1*10^9; % Hz
NP=W_c*273*1.38*(10)^-23; % watts/m^2
gamma_R = G_tx*G_rx*(c/(4*pi*f_R))^2;
lambda_i = 1;
lambda_s = 1;

% gamma_bar = (gamma_R*P_R)/NP;
mu = 0.002:0.0001:mu_max; % BS densities
L = 2000; % Length of road
d_safe = 5; % BS safety distance
h_bs = 8; % BS height


for i=1:length(mu)
    N = mu(i)*L;
    for j = 1:N-1
        d1(j) = (h_bs^2 + d_safe^2 + ((2*j+1)^2/(4*mu(i)^2)))^(-alpha/2);
    end
    lambda_ii1 = lambda_i./(gamma_R.*d1);

    rho_1 = ((1./(2.*mu(i))).^2 + h_bs^2 + d_safe^2)^(-alpha/2);

    lambda_rho1 = lambda_s./(gamma_R.*rho_1);

    %     fun1 = @(z)(exp(-NP*z).*(M_I(z, mu(i), L, gamma_R, lambda_i, d_safe(1), h_bs(1), alpha).*(1-M_S(z, mu(i), d_safe(1), h_bs(1), gamma_R, lambda_s, alpha)))./(z));
    R_avg1(i) = W_c.*(1-h_d.*mu(i).*V(1)).*(1/(log(2))).*closed(lambda_ii1, lambda_rho1);

    %     fun2 = @(z)(exp(-NP*z).*(M_I(z, mu(i), L, gamma_R, lambda_i, d_safe(2), h_bs(2), alpha).*(1-M_S(z, mu(i), d_safe(2), h_bs(2), gamma_R, lambda_s, alpha)))./(z));
    R_avg2(i) = W_c.*(1-h_d.*mu(i).*V(2)).*(1/(log(2))).*closed(lambda_ii1, lambda_rho1);

    %     fun3 = @(z)(exp(-NP*z).*(M_I(z, mu(i), L, gamma_R, lambda_i, d_safe(3), h_bs(3), alpha).*(1-M_S(z, mu(i), d_safe(3), h_bs(3), gamma_R, lambda_s, alpha)))./(z));
    R_avg3(i) = W_c.*(1-h_d.*mu(i).*V(3)).*(1/(log(2))).*closed(lambda_ii1, lambda_rho1);

    %     fun4 = @(z)(exp(-NP*z).*(M_I(z, mu(i), L, gamma_R, lambda_i, d_safe(4), h_bs(4), alpha).*(1-M_S(z, mu(i), d_safe(4), h_bs(4), gamma_R, lambda_s, alpha)))./(z));
    R_avg4(i) = W_c.*(1-h_d.*mu(i).*V(4)).*(1/(log(2))).*closed(lambda_ii1, lambda_rho1);
    
    r = 0;
    for k = 1:length(lambda_ii1)
        pd = makedist('Exponential', 'mu', 1/lambda_ii1(k));
        r = r + random(pd, 100000, 1); % Add all exponentials from interfering BSs
    end
    pd_s = makedist('Exponential', 'mu', 1/lambda_rho1); % Subtract recieved signal exp dist.
    z = random(pd_s, 100000, 1)./r;
    R_sim1(i) = W_c.*(1-h_d.*mu(i).*V(1)).*(1.6/(log(2))).*mean(log(1+z));
    R_sim2(i) = W_c.*(1-h_d.*mu(i).*V(2)).*(1.6/(log(2))).*mean(log(1+z));
    R_sim3(i) = W_c.*(1-h_d.*mu(i).*V(3)).*(1.6/(log(2))).*mean(log(1+z));
    R_sim4(i) = W_c.*(1-h_d.*mu(i).*V(4)).*(1.6/(log(2))).*mean(log(1+z));
end

% for j=1:length(mu)
%     if (W_c*R_avg1(j) < R_th)
%         break;
%     else
%       V_data1(j) = (1./(h_d.*mu(j))).*(1-R_th./(W_c.*R_avg1(j)));
%     end
% end
%
% for j1=1:length(mu)
%     if (W_c*R_avg2(j1) < R_th)
%         break;
%     else
%       V_data2(j1) = (1./(h_d.*mu(j1))).*(1-R_th./(W_c.*R_avg2(j1)));
%     end
% end
%
% for j2=1:length(mu)
%     if (W_c*R_avg3(j2) < R_th)
%         break;
%     else
%       V_data3(j2) = (1./(h_d.*mu(j2))).*(1-R_th./(W_c.*R_avg3(j2)));
%     end
% end
%
% tau = 0.006;
% epsilon = 0.01;
% sigma_LN = 1;
% mu_LN = 0;
% V_safe = (1/tau)*exp(sigma_LN*sqrt(2)*erfinv(2*epsilon -1) + mu_LN);
%
% for k=1:length(mu)
%     V_opt1(k) = min([V_safe V_max V_data1(k)]);
%     V_opt2(k) = min([V_safe V_max V_data2(k)]);
%     V_opt3(k) = min([V_safe V_max V_data3(k)]);
%
%     % Simulation
%     pd = makedist('Lognormal', 'mu', mu_LN, 'sigma', sigma_LN);
%     t = truncate(pd,0,Inf);
%     r = random(t, 100000, 1);
%     r_mean(k) = mean(1./r);
%     Qsim1(k) = V_opt1(k).*r_mean(k);
%     Qsim2(k) = V_opt2(k).*r_mean(k);
%     Qsim3(k) = V_opt3(k).*r_mean(k);
% end

figure(1)
% Y = exp((sigma_LN^2 - 2.*mu_LN)/2);
% Q1 = Y.*V_opt1;
% Q2 = Y.*V_opt2;
% Q3 = Y.*V_opt3;
semilogy(mu, R_avg1, 'k', 'LineWidth', 1.2)
hold on
semilogy(mu, R_avg2, 'b', 'LineWidth', 1.2)
hold on
semilogy(mu, R_avg3, 'r', 'LineWidth', 1.2)
hold on
semilogy(mu, R_avg4, 'm', 'LineWidth', 1.2)
hold on
semilogy(mu(1:10:end), R_sim1(1:10:end), 'ko', 'LineWidth', 1.2)
hold on
semilogy(mu(1:10:end), R_sim2(1:10:end), 'bo', 'LineWidth', 1.2)
hold on
semilogy(mu(1:10:end), R_sim3(1:10:end), 'ro', 'LineWidth', 1.2)
hold on
semilogy(mu(1:10:end), R_sim4(1:10:end), 'mo', 'LineWidth', 1.2)
hold off
xlabel('BS Density (\mu) [BSs/m]')
ylabel('Ergodic Capacity [bps]')
grid on;
legend('V = 0 m/s (Ana.)', 'V = 5 m/s (Ana.)', 'V = 15 m/s (Ana.)', 'V = 25 m/s (Ana.)', 'V = 0 m/s (Sim.)', 'V = 5 m/s (Sim.)', 'V = 15 m/s (Sim.)', 'V = 25 m/s (Sim.)');

function c = closed(lambda_i, lambda_s)
sum = 0;
for i = 1:length(lambda_i)
    W = 1;
    for j = 1:length(lambda_i)
        if j == i
            continue;
        end
        W = W*(lambda_i(j)/(lambda_i(j)-lambda_i(i)));
    end
    x(i) = log(lambda_i(i)/lambda_s) + atan2(lambda_i(i),0);
    sum = sum + W.*x(i).*(lambda_i(i)./(lambda_s - lambda_i(i)));
end
c = -(sum);
end

function h = hypo(z, mu, L, gamma_R, lambda_i, xs, zh, alpha)
cdf = 0;
N = mu*L;
for i = 1:N-1
    d(i) = (xs^2 + zh^2 + ((2*i+1)^2/(4*mu^2)))^(-alpha/2);
end

for i = 1:N-1
    W = 1;
    for j = 1:N-1
        if j == i
            continue;
        end
        W = W*((lambda_i/(gamma_R*d(j)))/((lambda_i/(gamma_R*d(j)))-(lambda_i/(gamma_R*d(i)))));
    end
    cdf = cdf + W.*((lambda_i/(gamma_R*d(i)))./(lambda_i.*z+(lambda_i/(gamma_R*d(i)))));
end
h = 1 - cdf;
end


function m = M_I(z, mu, L, gamma_R, lambda_i, xs, zh, alpha)
m = 1;
N = mu*L;
for i=1:N-1
    d = (xs^2 + zh^2 + ((2*i+1)^2/(4*mu^2)))^(-alpha/2);
    m = m.*(1./(1+z.*(((gamma_R.*d)./(lambda_i)))));
end
end

function ms = M_S(z, mu, xs, zh, gamma_R, lambda_s, alpha)
rho_0 = sqrt((1./(2.*mu)).^2 + zh^2 + xs^2);
ms = 1./(1+z.*(((gamma_R*rho_0.^(-alpha))./(lambda_s))));
end