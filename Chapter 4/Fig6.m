clc;
clear all;
close all;

mu_max = 0.01; % Max BS density
V_max = 30; % m/s
R_th =  1.2*10^8; % B/sec
alpha = 3;

mu_opt = fminbnd(@(mu)Vdata(mu, V_max, mu_max, alpha, R_th), 0.001, mu_max);

W_c = 40*10^6; % Hz
P_R = 1; % Watts
h_d = 15/(mu_max*V_max); % Normalized handover rate
% SINR parameters
% alpha = [2, 4]; % path loss exponent
G_tx = 1; % Gain at transmitter
G_rx = 1; % Gain at reciever
c = 3*10^8; % m/s
f_R = 2.1*10^9; % Hz
NP=W_c*273*1.38*(10)^-23; % watts/m^2
gamma_R = G_tx*G_rx*(c/(4*pi*f_R))^2;
lambda_i = 1;
lambda_s = 1;

% gamma_bar = (gamma_R*P_R)/NP;
% mu = 0.0009:0.0001:0.01; % BS densities
L = 2000; % Length of road
d_safe = 5; % BS safety distance
h_bs = 8; % BS height

fun = @(z)(exp(-NP*z).*(M_I(z, mu_opt, L, gamma_R, lambda_i, d_safe, h_bs, alpha).*(1-M_S(z, mu_opt, d_safe, h_bs, gamma_R, lambda_s, alpha)))./(z));
R_avg = W_c*integral(fun, 0, inf);
V_data = (1./(h_d.*mu_opt)).*(1-R_th./(W_c.*R_avg));

tau = [0.005, 0.006, 0.007];
epsilon = linspace(0.0001, 0.02);
sigma_LN = 1;
mu_LN = 0;

for i=1:length(epsilon)
    V_safe1(i) = (1/tau(1))*exp(sigma_LN*sqrt(2)*erfinv(2*epsilon(i) -1) + mu_LN);
    V_safe2(i) = (1/tau(2))*exp(sigma_LN*sqrt(2)*erfinv(2*epsilon(i) -1) + mu_LN);
    V_safe3(i) = (1/tau(3))*exp(sigma_LN*sqrt(2)*erfinv(2*epsilon(i) -1) + mu_LN);

    V_opt1(i) = min([V_safe1(i) V_max V_data]);
    V_opt2(i) = min([V_safe2(i) V_max V_data]);
    V_opt3(i) = min([V_safe3(i) V_max V_data]);


    pd = makedist('Lognormal', 'mu', mu_LN, 'sigma', sigma_LN);
    t = truncate(pd,0,Inf);
    r = random(t, 100000, 1);
    r_mean(i) = mean(1./r);
    Qsim1(i) = V_opt1(i).*r_mean(i);
    Qsim2(i) = V_opt2(i).*r_mean(i);
    Qsim3(i) = V_opt3(i).*r_mean(i);
end
Y = exp((sigma_LN^2 - 2.*mu_LN)/2);
Q1 = Y.*V_opt1;
Q2 = Y.*V_opt2;
Q3 = Y.*V_opt3;

figure(1)
plot(epsilon, Q1, 'k', 'LineWidth', 1.2)
hold on
plot(epsilon, Q2, 'b', 'LineWidth', 1.2)
hold on
plot(epsilon, Q3, 'r', 'LineWidth', 1.2)
hold on
plot(epsilon(1:10:end), Qsim1(1:10:end), 'ko', 'LineWidth', 1.2)
hold on
plot(epsilon(1:10:end), Qsim2(1:10:end), 'bo', 'LineWidth', 1.2)
hold on
plot(epsilon(1:10:end), Qsim3(1:10:end), 'ro', 'LineWidth', 1.2)
hold off;
xlabel('Crash Intensity Level (\epsilon)')
ylabel('Optimal Traffic Flow (Q)')
grid on;
legend('\tau = 0.005 s (Ana.)', '\tau = 0.006 s (Ana.)', '\tau = 0.007 s (Ana.)','\tau = 0.005 s (Sim.)', '\tau = 0.006 s (Sim.)', '\tau = 0.007 s (Sim.)');

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
    x(i) = log(lambda_i(i)/lambda_s) + 1i*atan2(lambda_i(i),0);
    sum = sum + W.*x(i).*(lambda_i(i)./(lambda_s - lambda_i(i)));
end
c = -real(sum);
end


function f = Vdata(mu, V_max, mu_max, alpha, R_th)
W_c = 40*10^6; % Hz
P_R = 1; % Watts
h_d = 3/(mu_max*V_max); % Normalized handover rate
% SINR parameters
% alpha = [2, 4]; % path loss exponent
G_tx = 1; % Gain at transmitter
G_rx = 1; % Gain at reciever
c = 3*10^8; % m/s
f_R = 2.1*10^9; % Hz
NP=W_c*273*1.38*(10)^-23; % watts/m^2
gamma_R = G_tx*G_rx*(c/(4*pi*f_R))^2;
lambda_i = 1;
lambda_s = 1;

% gamma_bar = (gamma_R*P_R)/NP;
% mu = 0.0009:0.0001:0.01; % BS densities
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

    R_avg1(i) = (1/(log(2))).*closed(lambda_ii1, lambda_rho1);
%     fun1 = @(z)(exp(-NP*z).*(M_I(z, mu(i), L, gamma_R, lambda_i, xs, zh, alpha).*(1-M_S(z, mu(i), xs, zh, gamma_R, lambda_s, alpha)))./(z));
%     R_avg1(i) = integral(fun1, 0, inf);

end

f = -max(0, (1./(h_d.*mu)).*(1-R_th./(W_c.*R_avg1)));
end

function m = M_I(z, mu, L, gamma_R, lambda_i, xs, zh, alpha)
m = 1;
N = mu*L;
for i=1:(N-1)/2
    d = (xs^2 + zh^2 + ((2*i+1)^2/(4*mu^2)))^(-alpha/2);
    m = m.*(1./(1+z.*(((gamma_R.*d)./(lambda_i)))));
end
end

function ms = M_S(z, mu, xs, zh, gamma_R, lambda_s, alpha)
rho_0 = sqrt((1./(2.*mu)).^2 + zh^2 + xs^2);
ms = 1./(1+z.*(((gamma_R*rho_0.^(-alpha))./(lambda_s))));
end