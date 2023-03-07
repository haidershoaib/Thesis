clc;
clear all;
close all;
V_max = 30; % m/s
% R_th =  5*10^7; % B/sec=
W_c = 40*10^6; % Hz
P_R = 1; % Watts
mu_max = 0.01; % Max BS density
h_d = 3/(mu_max*V_max); % Normalized handover rate

% SINR parameters
% alpha = [2, 4]; % path loss exponent
alpha = 2.5;
alpha1 = 3;
alpha2 = 3.5;
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
% mu = 0.001:0.0001:0.01; % BS densities
L = 2000; % Length of road
d_safe = 5; % BS safety distance
h_bs = 8; % BS height

R_th =  linspace(6*10^7,8*10^7);
% R_th = 6*10^7:2000000:8*10^7;
% R_th2 =  linspace(6*10^3,8*10^6);

for i = 1:length(R_th)
    tic
    mu_opt(i) = fminbnd(@(mu)Vdata(mu, R_th(i), alpha), 0.0015, 0.01);
    time(i) = toc;
end

for i = 1:length(R_th)
    mu_opt2(i) = fminbnd(@(mu)Vdata(mu, R_th(i), alpha1), 0.0015, 0.01);
end

for i = 1:length(R_th)
    mu_opt3(i) = fminbnd(@(mu)Vdata(mu, R_th(i), alpha2), 0.0015, 0.01);
end

mu = mu_opt;
mu1 = mu_opt2;
mu2 = mu_opt3;

% mu2 = mu_opt2;
for i=1:length(R_th)
    %     N = mu(i)*L;
    %     for j = 1:N-1
    %         d(j) = (h_bs^2 + d_safe^2 + ((2*j+1)^2/(4*mu(i)^2)))^(-alpha/2);
    %         d1(j) = (h_bs^2 + d_safe^2 + ((2*j+1)^2/(4*mu(i)^2)))^(-alpha1/2);
    %         d2(j) = (h_bs^2 + d_safe^2 + ((2*j+1)^2/(4*mu(i)^2)))^(-alpha2/2);
    %     end
    %     lambda = lambda_s./(gamma_R.*d);
    %
    %     rho_0 = ((1./(2.*mu(i))).^2 + h_bs^2 + d_safe^2)^(-alpha/2);
    %     rho_1 = ((1./(2.*mu(i))).^2 + h_bs^2 + d_safe^2)^(-alpha1/2);
    %     rho_2 = ((1./(2.*mu(i))).^2 + h_bs^2 + d_safe^2)^(-alpha2/2);
    %
    %     lambda2 = lambda_s./(Pj*gamma_R.*rho_0);
    %     lambda21 = lambda_s./(Pj*gamma_R.*rho_1);
    %     lambda22 = lambda_s./(Pj*gamma_R.*rho_2);

    fun = @(z)(exp(-NP*z).*(M_I(z, mu(i), L, Pj*gamma_R, lambda_i, d_safe, h_bs, alpha).*(1-M_S(z, mu(i), d_safe, h_bs, gamma_R, lambda_s, alpha)))./(z));
    R_avg(i) = integral(fun, 0, inf);

    fun1 = @(z)(exp(-NP*z).*(M_I(z, mu1(i), L, Pj*gamma_R, lambda_i, d_safe, h_bs, alpha1).*(1-M_S(z, mu1(i), d_safe, h_bs, gamma_R, lambda_s, alpha1)))./(z));
    R_avg1(i) = integral(fun1, 0, inf);
    
    fun2 = @(z)(exp(-NP*z).*(M_I(z, mu2(i), L, Pj*gamma_R, lambda_i, d_safe, h_bs, alpha2).*(1-M_S(z, mu2(i), d_safe, h_bs, gamma_R, lambda_s, alpha2)))./(z));
    R_avg2(i) = integral(fun1, 0, inf);

    %     fun2 = @(z)((1-hypo(z, lambda, lambda21))./(1+z));
    %     erg_cap2(i) = W_c.*(1/(log(2))).*integral(fun2, 0, inf);
end
V_data = (1./(h_d.*mu)).*(1-R_th./(W_c.*R_avg));
V_data1 = (1./(h_d.*mu1)).*(1-R_th./(W_c.*R_avg1));
V_data2 = (1./(h_d.*mu2)).*(1-R_th./(W_c.*R_avg2));

% for j=1:length(mu)
%     if (W_c*R_avg(j) < R_th)
%         break;
%     else
%       V_data(j) = (1./(h_d.*mu(j))).*(1-R_th./(W_c.*R_avg(j)));
%       V_data1(j) = (1./(h_d.*mu(j))).*(1-R_th./(W_c.*R_avg1(j)));
%     end
% end
tau = 0.006;
epsilon = 0.01;
sigma_LN = 1;
mu_LN = 0;
V_safe = (1/tau)*exp(sigma_LN*sqrt(2)*erfinv(2*epsilon -1) + mu_LN);

for k=1:length(mu)
    V_opt1(k) = min([V_safe V_max V_data(k)]);
    V_opt2(k) = min([V_safe V_max V_data1(k)]);
    V_opt3(k) = min([V_safe V_max V_data2(k)]);

    % Simulation
    pd = makedist('Lognormal', 'mu', mu_LN, 'sigma', sigma_LN);
    t = truncate(pd,0,Inf);
    r = random(t, 100000, 1);
    r_mean(k) = mean(1./r);
    Qsim1(k) = V_opt1(k).*r_mean(k);
    Qsim2(k) = V_opt2(k).*r_mean(k);
    Qsim3(k) = V_opt3(k).*r_mean(k);
end

figure(1)
subplot(2,1,1)
Y = exp((sigma_LN^2 - 2.*mu_LN)/2);
Q1 = Y.*V_opt1;
Q2 = Y.*V_opt2;
Q3 = Y.*V_opt3;
plot(R_th, Q1, 'k', 'LineWidth', 1.2);
hold on
plot(R_th, Q2, 'b', 'LineWidth', 1.2);
hold on
plot(R_th, Q3, 'r', 'LineWidth', 1.2);
hold on
plot(R_th(1:5:end), Qsim1(1:5:end), 'ko', 'LineWidth', 1.2);
hold on
plot(R_th(1:5:end), Qsim2(1:5:end), 'bo', 'LineWidth', 1.2);
hold on
plot(R_th(1:5:end), Qsim3(1:5:end), 'ro', 'LineWidth', 1.2);
grid on
xlabel('R_{th} [bps]');
ylabel('Optimal Traffic Flow (Q)');
legend('\alpha = 2.5 (Ana.)', '\alpha = 3 (Ana.)', '\alpha = 3.5 (Ana.)', '\alpha = 2.5 (Sim.)', '\alpha = 3 (Sim.)', '\alpha = 3.5 (Sim.)')

subplot(2,1,2)
plot(R_th, mu, '-ko', 'MarkerIndices',1:5:length(mu), 'LineWidth', 1.2);
hold on
plot(R_th, mu1, '-bo', 'MarkerIndices',1:5:length(mu1), 'LineWidth', 1.2);
hold on
plot(R_th, mu2, '-ro', 'MarkerIndices',1:5:length(mu2), 'LineWidth', 1.2);
grid on
legend('\alpha = 2.5', '\alpha = 3', '\alpha = 3.5')
xlabel('R_{th} [bps]');
ylabel('Optimal BS Density (\mu^*)');
avgT = mean(time)
function f = Vdata(mu, R_th, alpha)
V_max = 30; % m/s
% R_th =  5*10^7; % B/sec
W_c = 40*10^6; % Hz
P_R = 1; % Watts
mu_max = 0.01; % Max BS density
h_d = 3/(mu_max*V_max); % Normalized handover rate

% SINR parameters
% alpha = [2, 4]; % path loss exponent
% alpha = 3;
% alpha1 = 2;
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
% mu = 0.001:0.0001:0.01; % BS densities
L = 2000; % Length of road
d_safe = 5; % BS safety distance
h_bs = 8; % BS height
for i=1:length(mu)
    N = mu(i)*L;
    for j = 1:N-1
        d(j) = (h_bs^2 + d_safe^2 + ((2*j+1)^2/(4*mu(i)^2)))^(-alpha/2);
        %         d1(j) = (h_bs^2 + d_safe^2 + ((2*j+1)^2/(4*mu(i)^2)))^(-alpha1/2);
    end
    lambda = lambda_s./(gamma_R.*d);

    rho_0 = ((1./(2.*mu(i))).^2 + h_bs^2 + d_safe^2)^(-alpha/2);
    %     rho_1 = ((1./(2.*mu(i))).^2 + h_bs^2 + d_safe^2)^(-alpha1/2);
    lambda2 = lambda_s./(Pj*gamma_R.*rho_0);
    %     lambda21 = lambda_s./(Pj*gamma_R.*rho_1);

    fun = @(z)(exp(-NP*z).*(M_I(z, mu(i), L, Pj*gamma_R, lambda_i, d_safe, h_bs, alpha).*(1-M_S(z, mu(i), d_safe, h_bs, gamma_R, lambda_s, alpha)))./(z));
    R_avg(i) = integral(fun, 0, inf);

    %     fun1 = @(z)(exp(-NP*z).*(M_I(z, mu(i), L, Pj*gamma_R, lambda_i, d_safe, h_bs, alpha1).*(1-M_S(z, mu(i), d_safe, h_bs, gamma_R, lambda_s, alpha1)))./(z));
    %     R_avg1(i) = integral(fun1, 0, inf);


    %     fun2 = @(z)((1-hypo(z, lambda, lambda21))./(1+z));
    %     erg_cap2(i) = W_c.*(1/(log(2))).*integral(fun2, 0, inf);
end

for j=1:length(mu)
    if (W_c*R_avg(j) < R_th)
        break;
    else
        V_data(j) = (1./(h_d.*mu(j))).*(1-R_th./(W_c.*R_avg(j)));
        %       V_data1(j) = (1./(h_d.*mu(j))).*(1-R_th./(W_c.*R_avg1(j)));
    end
end
f = -V_data;

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