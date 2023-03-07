clc;
clear all;
close all;
v_max = 30; % m/s
R_thr = 1*10^9; % B/sec
W_c = 40*10^6; % Hz
P_R = 1; % Watts
%mu = 0.0005:0.000005:0.050; % BS intensity
mu_max = 1; % Max BS intensity
h_d = 1/(0.45*v_max); % Normalized handover rate

 
% SNR parameters
alpha = 3; % path loss exponent
G_tx = 1; % Gain at transmitter
G_rx = 1; % Gain at reciever
c = 3*10^8; % m/s
f_R = 2.1*10^9; % Hz
NP=W_c*273*1.38*(10)^-23; % watts/m^2
gamma_R = G_tx*G_rx*(c/(4*pi*f_R))^2;

gamma_bar = (gamma_R*P_R)/NP;

mu_data = (1/(2*gamma_bar^(1/alpha)))*nthroot((2^(R_thr/W_c)-1), alpha);
mu = linspace(mu_data, mu_max);
% Rate calculation
rho_0 = 1./(2*mu); % Average distance between BS and CAV
SNR = (gamma_R.*P_R.*(rho_0).^(-alpha))/NP;
R_m = W_c.*log2(1+SNR);
R_m_test = W_c.*log2(1+SNR);        
% Velocity with data rate
V_data = 1./(h_d.*mu).*((1 - (R_thr./R_m)));

%V_min calculation
% C1 parameters
tau_1 = 0.00015;
tau_2 = 0.0004;

lambda_bar1 = 1;
lambda_bar2 = 1;

% Typical crash probabilities are very small, gnerally between 0.1% and up
% to 5% 
epsilon_1 = 0.0015; % Crash probability for V_min below peak
epsilon_2 = 0.0015; % Crash probability for V_min above peak

v_safe_1 = -lambda_bar1*log(1-epsilon_1)/(tau_1); %Safe speed for V_min below peak
v_safe_2 = -lambda_bar2*log(1-epsilon_2)/(tau_2); %Safe speed for V_min above peak

v_min_1 = min([v_safe_1, v_max]); %V_min below peak
v_min_2 = min([v_safe_2, v_max]); %V_min above peak

% Plots for V_data vs mu for different crash probabilties
subplot(2,1,1)
plot(mu, v_min_1+mu*0, 'LineWidth', 1.2);
hold on;
plot(mu, V_data(1,:), 'LineWidth', 1.2)
hold off
%title('V_{data} vs. \mu')
xlabel('\mu (Number of BSs/m)')
ylabel('V_{data} (m/s)')
ylim([0 11])
xline(mu_data)
xline(0.45)
grid on;
legend('V_{min} = min({V_{max}, V_{safe}})', 'V_{data}');

subplot(2,1,2)
plot(mu, v_min_2+mu*0, 'LineWidth', 1.2);
hold on;
plot(mu, V_data(1,:), 'LineWidth', 1.2)
hold off
%title('V_{data} vs. \mu')
xlabel('\mu (Number of BSs/m)')
ylabel('V_{data} (m/s)')
% ylim([0 23])
xline(mu_data)
xline(0.45)
grid on;
legend('V_{min} = min({V_{max}, V_{safe}})', 'V_{data}');