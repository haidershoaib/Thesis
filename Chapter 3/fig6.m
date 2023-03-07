clc;
clear all;
close all;

% System Parameters
V_max = 30; % m/s
%V_max = [0.05 0.1 0.15]; % km/s
R_th = 1*10^9; % B/s
%R_thr = [1*10^9, 1*10^10, 1.1*10^10];
W_c = 40*10^6; % Hz
P_R = 1; % Watts
mu_max = 0.45; % Max BS intensity
h_d = 1./(mu_max.*V_max); % Normalized handover rate
%l = [0.001, 0.005, 0.01, 0.015]; % CAV safety distances (km)
l = 1;
epsilon = 0.0015; % Crash probability
lambda_bar = 1; % CAV average spacing
lambda = 1./lambda_bar;
tau = [0.0003, 0.0004, 0.0005]; % Processing time
 
% SNR parameters
alpha = 3; % path loss exponent
G_tx = 1; % Gain at transmitter
G_rx = 1; % Gain at reciever
c = 3*10^8; % km/s
f_R = 2.1*10^9; % Hz
NP=W_c*273*1.38*(10)^-23; % watts/m^2
gamma_R = G_tx*G_rx*(c/(4*pi*f_R))^2;
gamma_bar = (gamma_R.*P_R)/NP;

% Optimal BS intensity
%mu = (1/2).*(nthroot((2.^((sqrt(R_thr).*sqrt(R_thr+4.*W_c.*alpha)+R_thr)./(2.*W_c)))./(gamma_bar), alpha)); 
mu_data = (1/(2*gamma_bar^(1/alpha)))*nthroot((2^(R_th/W_c)-1), alpha);
mu = linspace(mu_data, mu_max); % BS intensity

% Vector allocation
SNR = zeros(1, length(mu));
R_m = zeros(1, length(mu));
V_data1 = zeros(1, length(mu));
V_data2 = zeros(1, length(mu));
V_data3 = zeros(1, length(mu));

V_safe1 = zeros(1, length(mu));
V_safe2 = zeros(1, length(mu));
V_safe3 = zeros(1, length(mu));

V_opt1 = zeros(1, length(mu));
V_opt2 = zeros(1, length(mu));
V_opt3 = zeros(1, length(mu));

% Safety Velocity
V_safe1 = -lambda_bar*log(1-epsilon)/(tau(1)); 
V_safe2 = -lambda_bar*log(1-epsilon)/(tau(2));
V_safe3 = -lambda_bar*log(1-epsilon)/(tau(3));

    for j = 1:length(mu)
        %SNR calculation
        SNR(j) = (gamma_R.*P_R.*(1/(2.*mu(j))).^(-alpha))./NP;
       
        % Rate calculation
        R_m(j) = W_c*log2(1 + SNR(j));
       
        % Velocity with data rate
        V_data1(j) = (1./(h_d(1).*mu(j))).*(1 - (R_th./R_m(j)));
        %V_data2(j) = (1./(h_d(2).*mu(j))).*(1 - (R_thr./R_m(j)));
        %V_data3(j) = (1./(h_d(3).*mu(j))).*(1 - (R_thr./R_m(j)));
        
        % Optimal Velocity
        V_opt1(j) = min([V_safe1 V_max V_data1(j)]);
        V_opt2(j) = min([V_safe2 V_max V_data1(j)]);
        V_opt3(j) = min([V_safe3 V_max V_data1(j)]);
        
        % Simulation
        pd = makedist('Exponential', 'mu', lambda_bar);
        t = truncate(pd,l,Inf);
        r = random(t, 10000, 1);
        Qsim1(j) = V_opt1(j).*mean(1./r);
        Qsim2(j) = V_opt2(j).*mean(1./r);
        Qsim3(j) = V_opt3(j).*mean(1./r);
    end

% Plots    
figure(1);
Y = exp(l.*lambda).*expint(l.*lambda);
Q1 = V_opt1.*lambda.*Y;
plot(mu, Q1, 'LineWidth', 1, 'LineWidth', 1.2);
hold on
% Y2 = exp(l(2).*lambda).*gammainc(lambda,l(2).*lambda,'upper');
Q2 = V_opt2.*lambda.*Y;
plot(mu, Q2, 'LineWidth', 1, 'LineWidth', 1.2);
hold on 
% Y3 = exp(l(3).*lambda).*gammainc(lambda,l(3).*lambda,'upper');
Q3 = V_opt3.*lambda.*Y;
plot(mu, Q3, 'LineWidth', 1, 'LineWidth', 1.2);
hold on;
plot(mu(1:4:end), Qsim1(1:4:end), 'o')
hold on;
plot(mu(1:4:end), Qsim2(1:4:end), 'o')
hold on;
plot(mu(1:4:end), Qsim3(1:4:end), 'o')
hold off;
grid on
xlabel('BS Density (\mu)');
ylabel('Traffic Flow (Q)');
%ylim([0 1.5*10^-4])
% title('Q vs. \mu with Varying Data Rate Thresholds');
legend('\tau = 0.00015 s (Ana.)', '\tau = 0.00020 s (Ana.)', '\tau = 0.00025 s (Ana.)');