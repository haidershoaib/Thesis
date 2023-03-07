clc;
clear all;
close all;

% System Parameters
V_max = 30; % m/s
%R_thr = 1*10^10; % B/s
% R_th = 1*10^9;

R_th = linspace(0.9*10^9, 1.1*10^9); 
W_c = 40*10^6; % Hz
P_R = 1; % Watts
mu_max = 0.45; % Max BS intensity
h_d = 1/(mu_max*V_max); % Normalized handover rate
l = 1;
%l = [0.001, 0.005, 0.01, 0.015]; % CAV safety distances (km)
epsilon = 0.0015; % Crash probability
lambda_bar = 1;
lambda = 1/lambda_bar;
%lambda = linspace(l, 5); % CAV average spacing
tau = 0.0002; % Processing time
 
% SNR parameters
alpha = [2, 3, 4]; % path loss exponent
G_tx = 1; % Gain at transmitter
G_rx = 1; % Gain at reciever
c = 3*10^8; % m/s
f_R = 2.1*10^9; % Hz
%NP=(10)^-23; % watts/m^2
NP=W_c*273*1.38*(10)^-23; % watts/m^2
gamma_R = G_tx*G_rx*(c/(4*pi*f_R))^2;
gamma_bar = (gamma_R.*P_R)/NP;

% Optimal BS intensity
% mu1 = (1/2).*(nthroot((2.^((sqrt(R_th).*sqrt(R_th+4.*W_c.*alpha(1))+R_th)./(2.*W_c)))./(gamma_bar), alpha(1))); 
% mu_hat2 = (1/2).*(nthroot((2.^((sqrt(R_th).*sqrt(R_th+4.*W_c.*alpha(2))+R_th)./(2.*W_c)))./(gamma_bar), alpha(2))); 
% mu_hat3 = (1/2).*(nthroot((2.^((sqrt(R_th).*sqrt(R_th+4.*W_c.*alpha(3))+R_th)./(2.*W_c)))./(gamma_bar), alpha(3)));
mu_hat1 = (1/(2*gamma_bar^(1/alpha(1)))).*(2.^((sqrt(R_th).*sqrt(R_th+4.*W_c.*alpha(1))+R_th)./(2.*W_c))).^(1/alpha(1)); 
mu_hat2 = (1/(2*gamma_bar^(1/alpha(2)))).*(2.^((sqrt(R_th).*sqrt(R_th+4.*W_c.*alpha(2))+R_th)./(2.*W_c))).^(1/alpha(2)); 
mu_hat3 = (1/(2*gamma_bar^(1/alpha(3)))).*(2.^((sqrt(R_th).*sqrt(R_th+4.*W_c.*alpha(3))+R_th)./(2.*W_c))).^(1/alpha(3));

mu_data1 = (1./(2.*gamma_bar.^(1./alpha(1))))*nthroot((2.^(R_th./W_c)-1), alpha(1));
mu_data2 = (1./(2.*gamma_bar.^(1./alpha(2))))*nthroot((2.^(R_th./W_c)-1), alpha(2));
mu_data3 = (1./(2.*gamma_bar.^(1./alpha(3))))*nthroot((2.^(R_th./W_c)-1), alpha(3));

mu_opt1 = zeros(1, length(R_th));
mu_opt2 = zeros(1, length(R_th));
mu_opt3 = zeros(1, length(R_th));

flag1 = zeros(1, length(R_th));
flag2 = zeros(1, length(R_th));
flag3 = zeros(1, length(R_th));
for k=1:length(R_th)
    if mu_hat1(k) < mu_data1(k)
       flag1(k) = 1;
    elseif mu_hat1(k) > mu_max
       mu_opt1(k) = mu_max;  
    else
        mu_opt1(k) = mu_hat1(k);
    end
    if mu_hat2(k) < mu_data2(k)
       flag2(k) = 1;
    elseif mu_hat2(k) > mu_max
       mu_opt2(k) = mu_max;
    else
        mu_opt2(k) = mu_hat2(k);
    end
    if mu_hat3(k) < mu_data3(k)
       flag3(k) = 1;
    elseif mu_hat3(k) > mu_max
       mu_opt3(k) = mu_max;
    else
         mu_opt3(k) = mu_hat3(k);
    end
end

% for i=1:length(R_th)
%     if mu_opt1(i) > mu_max
%        mu_opt1(i) = mu_max; 
%     end
%     if mu_opt2(i) > mu_max
%        mu_opt2(i) = mu_max; 
%     end
%     if mu_opt3(i) > mu_max
%        mu_opt3(i) = mu_max; 
%     end
% end


% mu4 = (1/2).*(nthroot((2.^((sqrt(R_thr).*sqrt(R_thr+4.*W_c.*alpha(4))+R_thr)./(2.*W_c)))./(gamma_bar), alpha(4)));
% mu = [mu1, mu2, mu3, mu4];
% Vector allocation
SNR1 = zeros(1, length(R_th));
SNR2 = zeros(1, length(R_th));
SNR3 = zeros(1, length(R_th));
% SNR4 = zeros(1, length(R_thr));

R_m1 = zeros(1, length(R_th));
R_m2 = zeros(1, length(R_th));
R_m3 = zeros(1, length(R_th));
% R_m4 = zeros(1, length(R_thr));

V_data1 = zeros(1, length(R_th));
V_data2 = zeros(1, length(R_th));
V_data3 = zeros(1, length(R_th));
% V_data4 = zeros(1, length(R_thr));

V_safe = zeros(1, length(R_th));

V_opt1 = zeros(1, length(R_th));
V_opt2 = zeros(1, length(R_th));
V_opt3 = zeros(1, length(R_th));
% V4 = zeros(1, length(R_thr));
for j = 1:length(R_th)
        V_safe(j) = -lambda_bar*log(1-epsilon)/(tau);
        if mu_opt1(j) > 0
            SNR1(j) = (gamma_R.*P_R.*((2.*mu_opt1(j))).^(alpha(1)))./NP;
            R_m1(j) = W_c*log2(1 + SNR1(j));
            V_data1(j) = (1/(h_d*mu_opt1(j)))*(1 - (R_th(j)./R_m1(j)));
            V_opt1(j) = min([V_safe(j) V_max max(V_data1(j), 0)]);
        else
            V_opt1(j) = 0;
        end
        
        if mu_opt2(j) > 0
            SNR2(j) = (gamma_R.*P_R.*((2.*mu_opt2(j))).^(alpha(2)))./NP;
            R_m2(j) = W_c*log2(1 + SNR2(j));
            V_data2(j) = (1/(h_d*mu_opt2(j)))*(1 - (R_th(j)./R_m2(j)));
            V_opt2(j) = min([V_safe(j) V_max max(V_data2(j), 0)]);
        else
            V_opt2(j) = 0;
        end
        
        if mu_opt3(j) > 0
            SNR3(j) = (gamma_R.*P_R.*((2.*mu_opt3(j))).^(alpha(3)))./NP;
            R_m3(j) = W_c*log2(1 + SNR3(j));
            V_data3(j) = (1/(h_d*mu_opt3(j)))*(1 - (R_th(j)./R_m3(j)));
            V_opt3(j) = min([V_safe(j) V_max max(V_data3(j), 0)]);
        else
            V_opt3(j) = 0;
        end
        % Simulation
          pd = makedist('Exponential', 'mu', lambda_bar);
          t = truncate(pd,l,Inf);
          r = random(t, 10000, 1);
          Qsim1(j) = V_opt1(j).*mean(1./r);
          Qsim2(j) = V_opt2(j).*mean(1./r);
          Qsim3(j) = V_opt3(j).*mean(1./r);
                  
end

% Plots    
% subplot(1,2,1)
figure(1)
% Y = exp(l.*lambda).*gamma(l.*lambda).*gammainc(lambda,l.*lambda,'upper');
Y = exp(l.*lambda).*expint(l.*lambda);
Q1 = V_opt1.*lambda.*Y;
plot(R_th, Q1, 'k', 'LineWidth', 1.2);
hold on
% Y2 = exp(l(2).*lambda).*gammainc(lambda,l(2).*lambda,'upper');
Q2 = V_opt2.*lambda.*Y;
plot(R_th, Q2, 'r', 'LineWidth', 1.2);
hold on 
% Y3 = exp(l(3).*lambda).*gammainc(lambda,l(3).*lambda,'upper');
Q3 = V_opt3.*lambda.*Y;
plot(R_th, Q3, 'b', 'LineWidth', 1.2);
hold on;
plot(R_th(1:4:end), Qsim1(1:4:end), 'ko')
hold on;
plot(R_th(1:4:end), Qsim2(1:4:end), 'ro')
hold on;
plot(R_th(1:4:end), Qsim3(1:4:end), 'bo')
hold off;
grid on
xlabel('Data Rate Threshold (R_{th}) [bps]');
ylabel('Traffic Flow (Q)');
% title('Q vs. \lambda with Varying Data Rate Thresholds');
legend('\alpha = 2 (Ana.)', '\alpha = 3 (Ana.)', '\alpha = 4 (Ana.)', '\alpha = 2 (Sim.)', '\alpha = 3 (Sim.)', '\alpha = 4 (Sim.)');
% figure(2)
% arr1 = ones(1, length(R_thr));
% r = V1./exprnd(arr1)/V_max;
% [f, x] = hist(r);
% bar(x, f)

% subplot(1,2,2)
% plot(R_th, mu_opt1, 'LineWidth', 1.2)
% hold on
% plot(R_th, mu_opt2, 'LineWidth', 1.2)
% hold on;
% plot(R_th, mu_opt3, 'LineWidth', 1.2)
% hold off;
% ylabel('Optimal BS Density (\mu^*)')
% xlabel('Data Rate Threshold (R_{th})')
% legend('\alpha = 2', '\alpha = 3', '\alpha = 4');
% grid on;

% figure(3)
% plot(R_th, V_opt1, 'LineWidth', 1.2)
% hold on
% plot(R_th, V_opt2, 'LineWidth', 1.2)
% hold on;
% plot(R_th, V_opt3, 'LineWidth', 1.2)
% hold off;
% ylabel('V_{data})')
% xlabel('Data Rate Threshold (R_{th})')
% legend('\alpha = 2', '\alpha = 3', '\alpha = 4');
% grid on;

% figure(4)
% plot(R_th, Qsim1)