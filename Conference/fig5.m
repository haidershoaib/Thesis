clc;
clear all;
close all;

% System Parameters
V_max = 30; % m/s
%R_thr = 1*10^10; % B/s
R_th = [1*10^9, 1.05*10^9 1.1*10^9]; % B/s
W_c = 40*10^6; % Hz
P_R = 1; % Watts
mu_max = 0.45; % Max BS intensity
h_d = 1/(mu_max*V_max); % Normalized handover rate
l = 1;
%l = [0.001, 0.005, 0.01, 0.015]; % CAV safety distances (km)
epsilon = 0.0015; % Crash probability
lambda_bar = linspace(l,10); % CAV average spacing
lambda = 1./lambda_bar;
tau = 0.0002; % Processing time
 
% SNR parameters
alpha = 3; % path loss exponent
G_tx = 1; % Gain at transmitter
G_rx = 1; % Gain at reciever
c = 3*10^8; % km/s
f_R = 2.1*10^9; % Hz
% NP=(10)^-23; % watts/m^2
NP=W_c*273*1.38*(10)^-23; % watts/m^2
gamma_R = G_tx*G_rx*(c/(4*pi*f_R))^2;
gamma_bar = (gamma_R.*P_R)/NP;

% Optimal BS intensity
mu_hat = (1/2).*(nthroot((2.^((sqrt(R_th).*sqrt(R_th+4.*W_c.*alpha)+R_th)./(2.*W_c)))./(gamma_bar), alpha)); 
mu_data = (1./(2.*gamma_bar.^(1./alpha)))*nthroot((2.^(R_th./W_c)-1), alpha);

mu_opt = zeros(1, length(mu_hat));
flag = zeros(1, length(mu_opt));
for k=1:length(mu_opt)
    if mu_hat(k) < mu_data(k)
       flag(k) = 1;
    elseif mu_hat(k) > mu_max
       mu_opt(k) = mu_max;  
    else
        mu_opt(k) = mu_hat(k);
    end
end

%mu = linspace(max([mu_data1 mu_data2 mu_data3]), mu_max); % BS intensity
% Vector allocation
SNR1 = zeros(1, length(lambda_bar));
SNR2 = zeros(1, length(lambda_bar));
SNR3 = zeros(1, length(lambda_bar));

R_m1 = zeros(1, length(lambda_bar));
R_m2 = zeros(1, length(lambda_bar));
R_m3 = zeros(1, length(lambda_bar));

V_data1 = zeros(1, length(lambda_bar));
V_data2 = zeros(1, length(lambda_bar));
V_data3 = zeros(1, length(lambda_bar));

V_safe = zeros(1, length(lambda_bar));

V_opt1 = zeros(1, length(lambda_bar));
V_opt2 = zeros(1, length(lambda_bar));
V_opt3 = zeros(1, length(lambda_bar));
    for j = 1:length(lambda_bar)
        V_safe(j) = -log(1-epsilon)/(lambda_bar(j)*tau); 
        if mu_opt(1) > 0
            SNR1(j) = (gamma_R.*P_R.*((2.*mu_opt(1))).^(alpha))./NP;
            R_m1(j) = W_c*log2(1 + SNR1(j));
            V_data1(j) = (1/(h_d*mu_opt(1)))*(1 - (R_th(1)./R_m1(j)));
            V_opt1(j) = min([V_safe(j) V_max V_data1(j)]);
        else
            V_opt1(j) = 0;
        end
        
        if mu_opt(2) > 0
            SNR2(j) = (gamma_R.*P_R.*((2.*mu_opt(2))).^(alpha))./NP;
            R_m2(j) = W_c*log2(1 + SNR2(j));
            V_data2(j) = (1/(h_d*mu_opt(2)))*(1 - (R_th(2)./R_m2(j)));
            V_opt2(j) = min([V_safe(j) V_max V_data2(j)]);
        else
            V_opt2(j) = 0;
        end
        
        if mu_opt(3) > 0
            SNR3(j) = (gamma_R.*P_R.*((2.*mu_opt(3))).^(alpha))./NP;
            R_m3(j) = W_c*log2(1 + SNR3(j));
            V_data3(j) = (1/(h_d*mu_opt(3)))*(1 - (R_th(3)./R_m3(j)));
            V_opt3(j) = min([V_safe(j) V_max V_data3(j)]);
        else
            V_opt3(j) = 0;
        end
        % Simulation
        pd = makedist('Exponential', 'mu', lambda(j));
        t = truncate(pd,l,Inf);
        r = random(t, 10000, 1);
        Qsim1(j) = V_opt1(j).*mean(1./r);
        Qsim2(j) = V_opt2(j).*mean(1./r);
        Qsim3(j) = V_opt3(j).*mean(1./r);
    end

% Plots    
figure(1);
Y = exp(l.*lambda_bar).*expint(l.*lambda_bar);
Q1 = V_opt1.*lambda_bar.*Y;
plot(lambda_bar, Q1, 'k', 'LineWidth', 1.2);
hold on
% Y2 = exp(l(2).*lambda).*gammainc(lambda,l(2).*lambda,'upper');
Q2 = V_opt2.*lambda_bar.*Y;
plot(lambda_bar, Q2, 'r', 'LineWidth', 1.2);
hold on 
% Y3 = exp(l(3).*lambda).*gammainc(lambda,l(3).*lambda,'upper');
Q3 = V_opt3.*lambda_bar.*Y;
plot(lambda_bar, Q3, 'b', 'LineWidth', 1.2);
hold on;
plot(lambda_bar(1:5:end), Qsim1(1:5:end), 'ko');
hold on;
plot(lambda_bar(1:5:end), Qsim2(1:5:end), 'ro');
hold on;
plot(lambda_bar(1:5:end), Qsim3(1:5:end), 'bo');
hold off;
grid on
xlabel('CAV Density (1/\lambda) [m^{-1}]');
ylabel('Traffic Flow (Q)');
% title('Q vs. 1/\lambda with Varying Data Rate Thresholds');
legend('R_{th} = 1 Gbps (Ana.)', 'R_{th} = 1.05 Gbps (Ana.)', 'R_{th} = 1.1 Gbps (Ana.)');

figure(2)
plot(lambda_bar, V_safe)
hold on;
plot(lambda_bar, V_data1+lambda_bar.*0, 'LineWidth', 1.2)
hold on;
plot(lambda_bar, V_data2+lambda_bar.*0, 'LineWidth', 1.2)
hold on;
plot(lambda_bar, V_data3+lambda_bar.*0, 'LineWidth', 1.2)
hold on;
plot(lambda_bar, V_max+lambda_bar.*0, 'LineWidth', 1.2)
hold off;
xlabel('CAV Density (1/\lambda) [m^{-1}]');
ylabel('Velocity in m/s');
legend('V_{safe}', 'V_{data} (R_{th} = 1 Gbps)', 'V_{data} (R_{th} = 1.05 Gbps)', 'V_{data} (R_{th} = 1.1 Gbps)', 'V_{max}');
grid on;

V = linspace(1, V_max);
R_mobil = R_m1.*(1 - h_d.*mu_opt(1).*V); 

figure(3)
plot(V, R_mobil)
grid on
ylabel('Data Rate [bps]');
xlabel('Speed [m/s]')