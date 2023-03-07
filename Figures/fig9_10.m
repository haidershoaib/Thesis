clear all; % This line clears all variables
close all;
anim=0; %binary variable equal to 1 if you want to see the animation, otherwise zero
rec=0; %binary variable equal to 1 if you want to record the animation
showFlow=0; %binary variable equal to 1 if you want to see the flow plot
gr=.1; %granularity of the simulation
TimeStep=40/gr; % T is the numebr of time steps (seconds) in the simulation
VehNum=20; % number of vehicles in the simulation
del=5; %initial distance between vehicles
% If initital distance is different for each vehicle, then the perimiter is
% calculated via summation
%animation properties
perim=15*del*1; 
rad=perim/2/pi;

x=linspace((VehNum-1)*del,0,VehNum); %the initial location of the N vehicles
V=ones(TimeStep,VehNum)*45*1000/3600; %initial speed of the vehicles
a=zeros(TimeStep,VehNum);

tau=0.1; 
%vMax=[45 50*ones(1,VehNum-1)]*1000/3600; %maximum speed
AcMax=5; %maximum acceleration m/s^2
DeMax=20; %maximum deceleration m/s^2

% Vector allocations for location of CAVs and dsitance from CAV to BS
x_cav = zeros(TimeStep-1, VehNum); 
y_cav = zeros(TimeStep-1, VehNum);
cavToBS = zeros(TimeStep-1, VehNum);

%*************************************************


% System Parameters
V_max = 30; %m/s
% R_th = 8*10^8; % B/s
P_R = 1; % Watts
mu_max = 0.45; % Max BS intensity
h_d = 1/(mu_max*V_max); % Normalized handover rate
l = 1;
tau_proccess = 0.4; %reaction of the drivers
for k=1:VehNum
   R_th(k) = 8*10^8; 
end

for k1=1:VehNum/2
   R_th(k1) = 8*10^8; 
end
% SNR parameters
alpha = 3; % path loss exponent
G_tx = 1; % Gain at transmitter
G_rx = 1; % Gain at reciever
c = 3*10^8; % m/s
f_R = 2.1*10^9; % Hz
W_c = 40*10^6;
NP=W_c*273*1.38*(10)^-23; % watts/m^2
gamma_R = G_tx*G_rx*(c/(4*pi*f_R))^2;
gamma_bar = (gamma_R.*P_R)/NP;


% Optimal BS intensity
mu_hat = (1/(2*gamma_bar^(1/alpha)))*(2^((sqrt(max(R_th))*sqrt(max(R_th)+4*W_c*alpha)+max(R_th))/(2*W_c)))^(1/alpha); 
mu_data = (1/(2*gamma_bar.^(1/alpha)))*nthroot((2^(max(R_th)/W_c)-1), alpha);
if mu_data > mu_hat
    return;
elseif mu_hat < mu_max
    BS_data = mu_data*perim; % Minimum
    BS_opt = mu_hat*perim; % Maximum
else
    BS_opt = mu_max*perim;
end

% Adding Equidistant BSs
angleTemp = 0;

% Uncomment below to get figure 10
% BS_opt = BS_data;
equalParts = 2*pi/ceil(BS_opt);
for j=1:ceil(BS_opt) 
   BSs(j) = angleTemp;
   angleTemp = angleTemp + equalParts;
end
mu_opt = length(BSs)/perim;
% BSs coordinates
x1 = rad*cos(BSs)+0; 
y1 = rad*sin(BSs)+0;

allDist = zeros(length(BSs));
arcDist = zeros(length(BSs));
%V_opt(1:TimeStep, 1:VehNum) = V_max;
%*************************************************
% The following "for" loops represents each time step of simulation using t

figure
for t=1:TimeStep-1
    for i=1:VehNum
        test(t, i) = x(i);
        %hMin=6+tau*V(c);
        hMin = V(t,i);
       if i==1
%             s(t, c)=abs(x(c)-x(VehNum)); % Distance between first vehicle and the last one
            V_opt(t, i) = V(t,i);
            h=mod(x(end),perim)-mod(x(1),perim);
            test1(t, i) = h;
            if h<0
                h=(h+perim)/tau_proccess;
                s(t, i) = h+perim;
            else
                s(t, i) = h;
            end
            if t>1
                SNR(t, i) = (exprnd(1)*gamma_R.*P_R.*(cavToBS(t-1, i)).^(-alpha))./NP;
                R_m(t, i) = W_c*log2(1 + SNR(t, i));
                V_data(t, i) = (1/(h_d*mu_opt))*(1 - (R_th(i)./R_m(t, i)));
                V_safe(t, i) = h;
                V_opt(t, i) = max(min([V_safe(t,i) V_data(t, i) V_max]), 0);
                h = V_opt(t, i);
            end
       else
            s(t, i) = x(i-1)-x(i); % Distance between vehicle and one behind it
            % At first time step, go at Vmax, no need to calculate distance
            % between CAV and BS yet
            if t ==1
                V_opt(t, i) = V_max;
            elseif t > 1
                % Calculate min of Vsafe, Vdata, Vmax
                SNR(t, i) = (exprnd(1)*gamma_R.*P_R.*(cavToBS(t-1, i)).^(-alpha))./NP;
                R_m(t, i) = W_c*log2(1 + SNR(t, i));
                V_data(t, i) = (1/(h_d*mu_opt))*(1 - (R_th(i)./R_m(t, i)));
                V_safe(t, i)= (s(t, i) - l)/tau_proccess;
                h = max(min([V_data(t, i) V_safe(t, i) V_max]), 0);
                V_opt(t, i) = h;
            end
        end
        % If V(t) > V_opt(t)
        if h<hMin
            V(t+1,i)=max(V(t,i)-tau*DeMax,0);
        % If V(t) < V_opt(t)     
        elseif h>hMin
            V(t+1,i)=min(V(t,i)+tau*AcMax,V_max);
        else
        % If V(t) = V_opt(t)
            V(t+1,i) = V(t,i);
            
        end
    end
    
%     if V(t+tau/gr,c)>V_max
%         V(t+tau/gr,c);
%     end
    
    length(V(t,:))
    x=x+V(t,:)*gr;
    X(t,:)=x;
    A(t,:)=(V(t+1,:)-V(t,:))/gr;
    xx1=mod(X(t,:),perim);
    teta1=xx1./rad;
    x_cav(t, :) = sin(teta1).*rad';
    y_cav(t, :) = cos(teta1).*rad';
    % Distance between CAV and BS
    for c1=1:VehNum
        cavToBS(t, c1) = dist(x_cav(t, c1), y_cav(t, c1), x1, y1, rad, length(BSs));
    end
end

%% Animation
if anim==1
    for t=1:TimeStep-1
        y=linspace(0,perim,1000);
        teta1=y/rad;
        for i=1:VehNum
            xx=mod(X(t,i),perim);
            teta=xx/rad;
%             x_cav1(t, c) = sin(teta)*rad';
%             y_cav1(t, c) = cos(teta)*rad';
%             % Distance between CAV and BS
%             cavToBS1(t, c) = dist(x_cav1(t, c), y_cav1(t, c), x1, y1, rad);
            h1=plot([sin(teta)*rad]',[cos(teta)*rad]','o');
            set(h1, 'markerfacecolor', get(h1, 'color'))
            hold on
        end
        plot(sin(teta1)*rad,cos(teta1)*rad,'k-');
        plot(x1, y1, 'o', 'LineWidth', 1) % Plot BSs on road
        axis equal
        axis off
        hold off
        if rec==1
            F(t) = getframe(gcf);
        else
            pause(.05)
        end
    end
    if rec==1
        writerObj = VideoWriter('PipeLowFastReaction.avi');
        writerObj.FrameRate = 10;
        % set the seconds per image
        % open the video writer
        open(writerObj);
        % write the frames to the video
        for i=1:length(F)
            % convert the image to a frame
            frame = F(i) ;
            writeVideo(writerObj, frame);
        end
        % close the writer object
        close(writerObj);
    end
else
    subplot(2,1,1)
    plot(X)
    ylabel('Distance')
    xlabel('Time')
    grid on;
    subplot(2,1,2)
    plot(V(1:TimeStep,:))
    ylabel('Speed')
    xlabel('Time')
    grid on;
    %axis equal;
%     subplot(3,1,3)
%     plot(A)
%     ylabel('Acceleration')
%     xlabel('Time')
%     grid on;
end

%% Traffic Flow
vAvg = mean(mean(V(1:end-1, :), 2));
% disp(mean(vAvg))
cav_density = VehNum/perim; 
flow = vAvg*cav_density;
disp(mean(flow))
if showFlow==1
    figure(2)
    plot(1:length(s), flow)
    xlabel('Time')
    ylabel('Traffic Flow')
end


%*************************************************
% Functions
function d = dist(vehLocX, vehLocY, baseLocX, baseLocY, rad, numBS)
    for p = 1:numBS
       allDist(p) = norm([vehLocX vehLocY] - [baseLocX(p) baseLocY(p)]);
       arcDist(p) = acos(1-allDist(p)^2/(2*rad));
    end
    d = min(allDist);
end