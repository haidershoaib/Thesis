clc;
clear all;
close all;

% System Parameters
V_max = 30; % m/s
%R_thr = 1*10^10; % B/s
R_th = [1*10^9, 1.05*10^9, 1.1*10^9];
W_c = 40*10^6; % Hz
P_R = 1; % Watts
mu_max = 0.45; % Max BS intensity
h_d = 1/(mu_max*V_max); % Normalized handover rate
%l = [0.001, 0.005, 0.01, 0.015]; % CAV safety distances (km)
l = 1;
epsilon = 0.0015; % Crash probability
lambda = 1; % CAV average spacing
tau = 0.00015; % Processing time
 
% SNR parameters
alpha = 3; % path loss exponent
G_tx = 1; % Gain at transmitter
G_rx = 1; % Gain at reciever
c = 3*10^8; % m/s
f_R = 2.1*10^9; % Hz
NP=W_c*273*1.38*(10)^-23; % watts/m^2
gamma_R = G_tx*G_rx*(c/(4*pi*f_R))^2;
gamma_bar = (gamma_R.*P_R)/NP;

% Optimal BS intensity
%mu = (1/2).*(nthroot((2.^((sqrt(R_thr).*sqrt(R_thr+4.*W_c.*alpha)+R_thr)./(2.*W_c)))./(gamma_bar), alpha)); 
mu_data1 = (1/(2*gamma_bar^(1/alpha)))*nthroot((2^(R_th(1)/W_c)-1), alpha);
mu_data2 = (1/(2*gamma_bar^(1/alpha)))*nthroot((2^(R_th(2)/W_c)-1), alpha);
mu_data3 = (1/(2*gamma_bar^(1/alpha)))*nthroot((2^(R_th(3)/W_c)-1), alpha);



if mu_data1 > mu_max
    flag1 = 1;
end
if mu_data2 > mu_max
    flag2 = 1;
end
if mu_data3 > mu_max
    flag3 = 1;
end

mu_opt = linspace(min([mu_data1 mu_data2 mu_data3]), mu_max); % BS intensity

% Vector allocation
SNR = zeros(1, length(mu_opt));
R_m = zeros(1, length(mu_opt));
V_data1 = zeros(1, length(mu_opt));
V_data2 = zeros(1, length(mu_opt));
V_data3 = zeros(1, length(mu_opt));

V_safe = zeros(1, length(mu_opt));

V_opt1 = zeros(1, length(mu_opt));
V_opt2 = zeros(1, length(mu_opt));
V_opt3 = zeros(1, length(mu_opt));
    for j = 1:length(mu_opt)
        %SNR calculation
        SNR(j) = (gamma_R.*P_R.*((2.*mu_opt(j))).^(alpha))./NP;
       
        % Rate calculation
        R_m(j) = W_c*log2(1 + SNR(j));
       
        % Velocity with data rate
        V_data1(j) = (1./(h_d.*mu_opt(j))).*(1 - (R_th(1)./R_m(j)));
        V_data2(j) = (1./(h_d.*mu_opt(j))).*(1 - (R_th(2)./R_m(j)));
        V_data3(j) = (1./(h_d.*mu_opt(j))).*(1 - (R_th(3)./R_m(j)));
      
        % Safety Velocity
        V_safe(j) = -log(1-epsilon)/(lambda*tau); %Safe speed for V_min below peak

        % Optimal Velocity
        V_opt1(j) = min([V_safe(j) V_max V_data1(j)]);
        V_opt2(j) = min([V_safe(j) V_max V_data2(j)]);
        V_opt3(j) = min([V_safe(j) V_max V_data3(j)]);
        
        % Simulation
        pd = makedist('Exponential', 'mu', lambda);
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
Q1 = max(Q1, 0);
plot(mu_opt, Q1, 'LineWidth', 1, 'LineWidth', 1.2);
hold on
Q2 = V_opt2.*lambda.*Y;
Q2 = max(Q2, 0);
plot(mu_opt, Q2, 'LineWidth', 1, 'LineWidth', 1.2);
hold on 
Q3 = V_opt3.*lambda.*Y;
Q3 = max(Q3, 0);
plot(mu_opt, Q3, 'LineWidth', 1, 'LineWidth', 1.2);
hold on;
Qsim1 = max(Qsim1, 0);
Qsim2 = max(Qsim2, 0);
Qsim3 = max(Qsim3, 0);
plot(mu_opt(2:4:end), Qsim1(2:4:end), 'o')
hold on;
plot(mu_opt(21:4:end), Qsim2(21:4:end), 'o')
hold on;
plot(mu_opt(48:4:end), Qsim3(48:4:end), 'o')
hold off;
grid on
xlabel('BS Density (\mu)');
ylabel('Traffic Flow (Q)');
% ylim([0 1.5*10^-4])
% title('Q vs. \mu with Varying Data Rate Thresholds');
% legend('R_{th} = 1 Gbps (Ana.)', 'R_{th} = 1.05 Gbps (Ana.)', 'R_{th} = 1.1 Gbps (Ana.)');