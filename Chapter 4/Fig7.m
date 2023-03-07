clc;
clear all;
close all;

V_max = 30; % m/s
R_th =  5.9*10^7; % B/sec
W_c = 40*10^6; % Hz
P_R = 1; % Watts
mu_max = 0.01; % Max BS density
h_d = 3/(mu_max*V_max); % Normalized handover rate

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
mu = 0.001:0.0001:0.01; % BS densities 
L = 2000; % Length of road
d_safe = [10, 70, 100]; % BS safety distance
h_bs = [50, 50, 50]; % BS height


for i=1:length(mu)
    N = mu(i)*L;
    for j = 1:N-1
        d1(j) = (h_bs(1)^2 + d_safe(1)^2 + ((2*j+1)^2/(4*mu(i)^2)))^(-alpha/2);
        d2(j) = (h_bs(2)^2 + d_safe(2)^2 + ((2*j+1)^2/(4*mu(i)^2)))^(-alpha/2);
        d3(j) = (h_bs(3)^2 + d_safe(3)^2 + ((2*j+1)^2/(4*mu(i)^2)))^(-alpha/2);
    end
    lambda_ii1 = lambda_s./(gamma_R.*d1);
    lambda_ii2 = lambda_s./(gamma_R.*d2);
    lambda_ii3 = lambda_s./(gamma_R.*d3);

    rho_1 = ((1./(2.*mu(i))).^2 + h_bs(1)^2 + d_safe(1)^2)^(-alpha/2);
    rho_2 = ((1./(2.*mu(i))).^2 + h_bs(2)^2 + d_safe(2)^2)^(-alpha/2);
    rho_3 = ((1./(2.*mu(i))).^2 + h_bs(3)^2 + d_safe(3)^2)^(-alpha/2);
    lambda_rho1 = lambda_s./(gamma_R.*rho_1);
    lambda_rho2 = lambda_s./(gamma_R.*rho_2);
    lambda_rho3 = lambda_s./(gamma_R.*rho_3);
    
    fun1 = @(z)((1-exponential(z, lambda_rho1))./(1+z));
    R_avg1_noint(i) = (1/(log(2))).*integral(fun1, 0, inf);
    R_avg1(i) = (1/(log(2))).*closed(lambda_ii1, lambda_rho1);
%     R_avg1(i) = (1/(log(2))).*integral(fun1, 0, inf);

    fun2 = @(z)((1-exponential(z, lambda_rho2))./(1+z));
    R_avg2_noint(i) = (1/(log(2))).*integral(fun2, 0, inf);
    R_avg2(i) = (1/(log(2))).*closed(lambda_ii2, lambda_rho2);
%     R_avg2(i) = (1/(log(2))).*integral(fun2, 0, inf);

    fun3 = @(z)((1-exponential(z, lambda_rho3))./(1+z));
    R_avg3_noint(i) = (1/(log(2))).*integral(fun3, 0, inf);
    R_avg3(i) = (1/(log(2))).*closed(lambda_ii3, lambda_rho3);
%     R_avg3(i) = (1/(log(2))).*integral(fun3, 0, inf);


%     fun1 = @(z)(exp(-NP*z).*(M_I(z, mu(i), L, gamma_R, lambda_i, d_safe(1), h_bs(1), alpha).*(1-M_S(z, mu(i), d_safe(1), h_bs(1), gamma_R, lambda_s, alpha)))./(z));
%     R_avg1(i) = integral(fun1, 0, inf);
% 
%     fun2 = @(z)(exp(-NP*z).*(M_I(z, mu(i), L, gamma_R, lambda_i, d_safe(2), h_bs(2), alpha).*(1-M_S(z, mu(i), d_safe(2), h_bs(2), gamma_R, lambda_s, alpha)))./(z));
%     R_avg2(i) = integral(fun2, 0, inf);
% 
%     fun3 = @(z)(exp(-NP*z).*(M_I(z, mu(i), L, gamma_R, lambda_i, d_safe(3), h_bs(3), alpha).*(1-M_S(z, mu(i), d_safe(3), h_bs(3), gamma_R, lambda_s, alpha)))./(z));
%     R_avg3(i) = integral(fun3, 0, inf);
end

for j0=1:length(mu)
    if (W_c*R_avg1(j0) < R_th)
        break;
    else
      V_data1(j0) = (1./(h_d.*mu(j0))).*(1-R_th./(W_c.*R_avg1(j0)));
    end
end

for j1=1:length(mu)
    if (W_c*R_avg2(j1) < R_th)
        break;
    else
      V_data2(j1) = (1./(h_d.*mu(j1))).*(1-R_th./(W_c.*R_avg2(j1)));
    end
end

for j2=1:length(mu)
    if (W_c*R_avg3(j2) < R_th)
        break;
    else
      V_data3(j2) = (1./(h_d.*mu(j2))).*(1-R_th./(W_c.*R_avg3(j2)));
    end
end
SINR_th = 2.^(R_th./W_c)-1; % SINR threshold

SINR_no_int1 = max(0, ((1/lambda_s)*P_R.*gamma_R.*rho_1)./(NP));
SINR_no_int1 = max(SINR_th, SINR_no_int1);

SINR_no_int2 = max(0, ((1/lambda_s)*P_R.*gamma_R.*rho_2)./(NP));
SINR_no_int2 = max(SINR_th, SINR_no_int2);

SINR_no_int3 = max(0, ((1/lambda_s)*P_R.*gamma_R.*rho_3)./(NP));
SINR_no_int3 = max(SINR_th, SINR_no_int3);

% Calculate V_data
R_m_no_int1 = W_c.*log2(1+SINR_no_int1);
V_data_no_int1 = max(0, 1./(h_d.*mu).*((1 - (R_th./R_m_no_int1))));

R_m_no_int2 = W_c.*log2(1+SINR_no_int2);
V_data_no_int2 = max(0, 1./(h_d.*mu).*((1 - (R_th./R_m_no_int2))));

R_m_no_int3 = W_c.*log2(1+SINR_no_int3);
V_data_no_int3 = max(0, 1./(h_d.*mu).*((1 - (R_th./R_m_no_int3))));


tau = 0.006;
epsilon = 0.01;
sigma_LN = 1;
mu_LN = 0;
V_safe = (1/tau)*exp(sigma_LN*sqrt(2)*erfinv(2*epsilon -1) + mu_LN);

for k=1:length(mu)
    V_opt1(k) = min([V_safe V_max V_data1(k)]);
    V_opt2(k) = min([V_safe V_max V_data2(k)]);
    V_opt3(k) = min([V_safe V_max V_data3(k)]);

    V_opt1_no_int(k) = min([V_safe V_max V_data_no_int1(k)]);
    V_opt2_no_int(k) = min([V_safe V_max V_data_no_int2(k)]);
    V_opt3_no_int(k) = min([V_safe V_max V_data_no_int3(k)]);

    % Simulation
    pd = makedist('Lognormal', 'mu', mu_LN, 'sigma', sigma_LN);
    t = truncate(pd,0,Inf);
    r = random(t, 100000, 1);
    r_mean(k) = mean(1./r);
    Qsim1(k) = V_opt1(k).*r_mean(k);
    Qsim2(k) = V_opt2(k).*r_mean(k);
    Qsim3(k) = V_opt3(k).*r_mean(k);

    Qsim1_no_int(k) = V_opt1_no_int(k).*r_mean(k);
    Qsim2_no_int(k) = V_opt2_no_int(k).*r_mean(k);
    Qsim3_no_int(k) = V_opt3_no_int(k).*r_mean(k);
end

figure(1)
Y = exp((sigma_LN^2 - 2.*mu_LN)/2);
Q1 = Y.*V_opt1;
Q2 = Y.*V_opt2;
Q3 = Y.*V_opt3;

Q1_no_int = Y.*V_opt1_no_int;
Q2_no_int = Y.*V_opt2_no_int;
Q3_no_int = Y.*V_opt3_no_int;

plot(mu, Q1, 'k', 'LineWidth', 1.2)
hold on
plot(mu, Q2, 'b', 'LineWidth', 1.2)
hold on
plot(mu, Q3, 'r', 'LineWidth', 1.2)
hold on
plot(mu, Q1_no_int, 'k--', 'LineWidth', 1.2)
hold on
plot(mu, Q2_no_int, 'b--', 'LineWidth', 1.2)
hold on
plot(mu, Q3_no_int, 'r--', 'LineWidth', 1.2)
hold on

plot(mu(1:4:end), Qsim1(1:4:end), 'ko', 'LineWidth', 1.2)
hold on
plot(mu(1:4:end), Qsim2(1:4:end), 'bo', 'LineWidth', 1.2)
hold on
plot(mu(1:4:end), Qsim3(1:4:end), 'ro', 'LineWidth', 1.2)
hold on


plot(mu(1:4:end), Qsim1_no_int(1:4:end), 'ko', 'LineWidth', 1.2)
hold on
plot(mu(1:4:end), Qsim2_no_int(1:4:end), 'bo', 'LineWidth', 1.2)
hold on
plot(mu(1:4:end), Qsim3_no_int(1:4:end), 'ro', 'LineWidth', 1.2)
hold off
xlabel('BS Density (\mu) [BSs/m]')
ylabel('Optimal Traffic Flow (Q)')
grid on;
legend('d_{safe} = 10 m (with interference)', 'd_{safe} = 70 m (with interference)', 'd_{safe} = 100 m (with interference)','d_{safe} = 10 m (w/o interference)', 'd_{safe} = 70 m (w/o interference)', 'd_{safe} = 100 m (w/o interference)');

% figure(2)
% Y = exp((sigma_LN^2 - 2.*mu_LN)/2);
% Q1 = Y.*V_opt1;
% Q2 = Y.*V_opt2;
% Q3 = Y.*V_opt3;
% plot(mu, Q1, 'k', 'LineWidth', 1.2)
% hold on
% plot(mu, Q2, 'b', 'LineWidth', 1.2)
% hold on
% plot(mu, Q3, 'r', 'LineWidth', 1.2)
% hold on
% plot(mu(1:4:end), Qsim1(1:4:end), 'ko', 'LineWidth', 1.2)
% hold on
% plot(mu(1:4:end), Qsim2(1:4:end), 'bo', 'LineWidth', 1.2)
% hold on
% plot(mu(1:4:end), Qsim3(1:4:end), 'ro', 'LineWidth', 1.2)
% hold off
% xlabel('BS Density (\mu) [BSs/m]')
% ylabel('Traffic Flow (Q)')
% grid on;
% legend('d_{safe}, h_{BS} = 1 m (Ana.)', 'd_{safe}, h_{BS} = 50 m (Ana.)', 'd_{safe}, h_{BS} = 70 m (Ana.)','d_{safe}, h_{BS} = 1 m (Sim.)', 'd_{safe}, h_{BS} = 50 m (Sim.)', 'd_{safe}, h_{BS} = 70 m (Sim.)');


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

function f = exponential(x, lambda)
    f = 1 - exp(-lambda*x);
end

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