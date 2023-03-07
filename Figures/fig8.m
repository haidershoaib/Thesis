clc
clear all; % This line clears all variables
close all;
anim=1; %binary variable equal to 1 if you want to see the animal, otherwise zero
rec=0; %binary variable equal to 1 if you want to record the animation
gr=.1; %granularity of the simulation
TimeStep=20/gr; % T is the numebr of time steps (seconds) in the simulation
VehNum=10; % number of vehicles in the simulation
del=10; %initial distance between vehicles

%animation properties
perim=VehNum*del*1.1;
rad=perim/2/pi;

x=linspace((VehNum-1)*del,0,VehNum); %the initial location of the 5 vehicles
V=ones(TimeStep,VehNum)*45*1000/3600; %initial speed of the vehicles
a=zeros(TimeStep,VehNum);

tau=.1; %reaction of the drivers
vMax=[45 50*ones(1,VehNum-1)]*1000/3600; %maximum speed
AcMax=5; %maximum acceleration m/s^2
DeMax=20; %maximum deceleration m/s^2

%*************************************************
% The following "for" loops represents each time step of simulation using t

figure
for t=1:TimeStep-1
    for c=1:VehNum
        hMin=6+tau*V(c);
        if c==1
            h=mod(x(end),perim)-mod(x(1),perim);
            if h<0
                h=h+perim;
            end
        else
            h=x(c-1)-x(c);
        end
        
        if h<hMin
            V(t+tau/gr,c)=max(V(t,c)-tau*DeMax,0);
        elseif h>hMin
            V(t+tau/gr,c)=min(V(t,c)+tau*AcMax,vMax(c));
            
        end
    end
    
    if V(t+tau/gr,c)>vMax
        V(t+tau/gr,c)
    end
    
    
    x=x+V(t,:)*gr;
    X(t,:)=x;
    A(t,:)=(V(t+1,:)-V(t,:))/gr;
end



if anim==1
    for t=1:TimeStep-1
        y=linspace(0,perim,1000);
        teta1=y/rad;
        for c=1:VehNum
            xx=mod(X(t,c),perim);
            teta=xx/rad;
            h1=plot([sin(teta)*rad]',[cos(teta)*rad]','o');
            set(h1, 'markerfacecolor', get(h1, 'color'))
            hold on
        end
        plot(sin(teta1)*rad,cos(teta1)*rad,'k-');
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
    grid on
    subplot(2,1,2)
    plot(V(1:TimeStep,:))
    ylabel('Speed')
    xlabel('Time')
    grid on
%     subplot(3,1,3)
%     plot(A)
%     ylabel('Acceleration')
%     xlabel('Time')
end