clc;
clear all;
close all;

mu_max = 0.02;
V_max = 30; % m/s

W_c = 40*10^6; % Hz
R_th = [1*10^8, 1*10^8, 1*10^8];
% gamma_th = 2.^(R_th./W_c)-1;
lambda = 1;
V = [5, 15, 25]; % m/s

P_R = 1;
alpha = 3;
G_tx = 1; % Gain at transmitter
G_rx = 1; % Gain at reciever
c = 3*10^8; % m/s
f_R = 2.1*10^9; % Hz
NP=W_c*273*1.38*(10)^-23; % watts/m^2
gamma_R = G_tx*G_rx*(c/(4*pi*f_R))^2;
x = 5; % BS safety distance
z = 8; % BS height
L = 2000; % Length of road
mu = 0.001:0.001:mu_max; % BS densities
% mu = linspace(0.1, mu_max);
for i = 1:length(mu)
    N = mu(i)*L;
    d_ij = (x^2 + z^2 + ((1)^2/(4*mu(i)^2)))^(1/2);

    for k = 1:N-1
        d(k) = (x^2 + z^2 + ((2*k+1)^2/(4*mu(i)^2)))^(1/2);
        lambda_i(k) = lambda/(P_R*gamma_R*d(k)^(-alpha));
    end
    if i==1
        test = lambda_i;
        test2 = lambda/(gamma_R*P_R*d_ij^(-alpha));
    end
    H_c = (mu(i).*V)./(mu_max.*V_max);
    gamma_th = 2.^(R_th(1)./(W_c))-1;
    gamma_th1 = 2.^(R_th(1)./(W_c*(1-H_c(1))))-1;
    gamma_th2 = 2.^(R_th(2)./(W_c*(1-H_c(2))))-1;
    gamma_th3 = 2.^(R_th(3)./(W_c*(1-H_c(3))))-1;

    lambda_s = lambda/(gamma_R*P_R*d_ij^(-alpha));

    rprob(i) = hypo(gamma_th, lambda_i, lambda_s);
    r_sim(i) = hypo_sim(gamma_th, lambda_i, lambda_s);

    r_prob1(i) = hypo(gamma_th1, lambda_i, lambda_s);
    r_sim1(i) = hypo_sim(gamma_th1, lambda_i, lambda_s);

    r_prob2(i) = hypo(gamma_th2, lambda_i, lambda_s);
    r_sim2(i) = hypo_sim(gamma_th2, lambda_i, lambda_s);

    r_prob3(i) = hypo(gamma_th3, lambda_i, lambda_s);
    r_sim3(i) = hypo_sim(gamma_th3, lambda_i, lambda_s);
end

figure(1)
plot(mu, rprob, '-ko', 'MarkerIndices',1:3:length(mu), 'LineWidth', 1.2)
hold on
plot(mu, r_prob1, '-bo', 'MarkerIndices',1:3:length(mu), 'LineWidth', 1.2)
hold on
plot(mu, r_prob2, '-ro', 'MarkerIndices', 1:3:length(mu), 'LineWidth', 1.2)
hold on
plot(mu, r_prob3, '-mo', 'MarkerIndices',1:3:length(mu), 'LineWidth', 1.2)
% hold on
% plot(mu(1:4:end), r_sim(1:4:end), 'ko', 'LineWidth', 1.2)
% hold on
% plot(mu(1:4:end), r_sim1(1:4:end), 'bo', 'LineWidth', 1.2)
% hold on
% plot(mu(1:4:end), r_sim2(1:4:end), 'ro', 'LineWidth', 1.2)
% hold on
% plot(mu(1:4:end), r_sim3(1:4:end), 'mo', 'LineWidth', 1.2)
hold off
xlabel('BS Density (\mu) [BS/m]')
ylabel('P_{out}')
grid on
legend('V = 0 m/s (Ana.)', 'V = 5 m/s (Ana.)', 'V = 15 m/s (Ana.)', 'V = 25 m/s (Ana.)', 'V = 0 m/s (Sim.)', 'V = 5 m/s (Sim.)', 'V = 15 m/s (Sim.)', 'V = 25 m/s (Sim.)')

% figure(2)
% plot(mu, r_test)

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

function r_simu = hypo_sim(z, lambda_i, lambda_s)
r = 0;
for i = 1:length(lambda_i)
    pd = makedist('Exponential', 'mu', 1/lambda_i(i));
    r = r + random(pd, 1000, 1); % Add all exponentials from interfering BSs
end
pd_s = makedist('Exponential', 'mu', 1/lambda_s);
rz = random(pd_s, 1000, 1)./r;

[f, x] = hist(rz, 100);
y = f./sum(f);
cdf_sim = cumsum(y);
r_simu = interp1(x, cdf_sim, z);
end