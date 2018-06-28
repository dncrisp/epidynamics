%% Simulations for Fig_S10
% Marisa Saggio
%
% Class SN/FLC, with three different paths giving different trends for the
% frequency behavior
% The model used is described in Saggio et al. 2017, JMN

%% 
clear all;

%------- SETTINGS - INTEGRATION
options = odeset('MaxStep', 0.1);   % maximal time integration step
dt=0.1;
tspan=[0:0.5:15000];                  % time span to integrate over
x=zeros(length(tspan),3,4);           % initial conditions (must be a column)

%------- SETTINGS MODEL
b = 1.0;                            % focus
R = 0.4;                            % radius in unfolding
c = 0.001;                          % velocity of slow variable
dstar =0.3 ;                        % threshold for slow variable inversion
N=1;                                % if N=1, upper branch of equilibria is resting state

%% INTEGRATION

FigHandle = figure('Position', [100, 100, 800, 450]);

for i=1:3
    
    % sets path parameters. p1 is called A in the paper, p2 is called B.
    switch i
        case 1
            p1=[0.2997,-0.1272,-0.9455]*0.4;p2=[0.7124,0.1464,-0.6863]*0.4;c = 0.001;  % constant frequency                 
            range_axis=[1700,3400,-1.2,0.7];
        case 2
            c=0.002; p1=[-0.2578,-0.2489,-0.9336]*0.4;  p2=[0.3591,0.08283,-0.1555]*0.4; % increasing frequency
            range_axis=[2400,4600,-1.2,0.7];
        case 3
            p1=[0.7046,0.05218,-0.7077]*0.4; p2=[0.5467,0.09841,-0.8315]*0.4; c=0.0003; % decreasing frequency
            range_axis=[2300,4000,-1,0.8];       
    end
    
    % computes path vectors E,F
    pin=[p1,p2];
    E=p1/R;
    F=cross(cross(p1,p2),p1);
    F=F/norm(F);
    
    % integrates equations
    [t,x(:,:,i)]=ode45('Cod3_Model',tspan, x0, options, b, c, R, dstar,E,F,N); % evaluation of stochastic ODE

    % ------- PLOT TIMESERIES
   
    subplot(3,1,i)
    str=sprintf('A = %.4f,%.4f,%.4f;  B = %.4f,%.4f,%.4f, c = %.5f, dstar = %.2f', pin(1), pin(2), pin(3), pin(4), pin(5), pin(6), c, dstar);
    title(str)
    newblue=[0, 113/255, 188/255];
    newred=[216/255, 82/255, 24/255];
    
        plot(tspan, x(:,1,i),'Color',newblue);
        hold on
        plot(tspan, x(:,3,i),'Color',newred);
        hold off
        title(str)
        %xlabel('time [s]')
        legend('x','z')
        axis(range_axis)
end

set(gcf, 'Color', 'w')

%% PLOT PATHS ON MAP 
%subplot(4,2,1)
openfig('../../0_COD3MODEL/Figures/Ampl_Equazim.fig')
hold on
topolegend=[1,90,0];
axesm('MapProjection','eqaazim','Grid','on')
for i=1:3
        switch i
       case 1
            p1=[0.2997,-0.1272,-0.9455]*0.4;p2=[0.7124,0.1464,-0.6863]*0.4;c = 0.001;  
        case 2
            c=0.002; p1=[-0.2578,-0.2489,-0.9336]*0.4;  p2=[0.3591,0.08283,-0.1555]*0.4;
        case 3
            p1=[0.7046,0.05218,-0.7077]*0.4; p2=[0.5467,0.09841,-0.8315]*0.4; c=0.0003; 
        end
        
        pin=[p1,p2];
        E=p1/R;
        F=cross(cross(p1,p2),p1);
        F=F/norm(F);
        
    % path followed on the map - cartesian coordinates for parameter space
        mu2=R*(E(1)*cos(x(:,3,i))+F(1)*sin(x(:,3,i)));
        mu1=-R*(E(2)*cos(x(:,3,i))+F(2)*sin(x(:,3,i)));
        nu=R*(E(3)*cos(x(:,3,i))+F(3)*sin(x(:,3,i)));
    % path followed on the map - spherical coordinates for parameter space
    ph = atan2(-mu1,mu2);
    th = acos(nu/R);
    % transform to latitude and longitude
    path_lat=real((pi/2-th)*180/pi);
    path_long=real(ph*180/pi);
    %plot
    geoshow(path_lat,path_long,'Color','w','LineWidth',2)
  axis([-.47,0.27,-1.57,-0.35])
end
set(gcf, 'Color', 'w')


openfig('../../0_COD3MODEL/Figures/Freq_Equazim.fig')
hold on
topolegend=[1,90,0];
axesm('MapProjection','eqaazim','Grid','on')
for i=1:3
        switch i
       case 1
            p1=[0.2997,-0.1272,-0.9455]*0.4;p2=[0.7124,0.1464,-0.6863]*0.4;c = 0.001; 
        case 2
            c=0.002; p1=[-0.2578,-0.2489,-0.9336]*0.4;  p2=[0.3591,0.08283,-0.1555]*0.4;
        case 3
            p1=[0.7046,0.05218,-0.7077]*0.4; p2=[0.5467,0.09841,-0.8315]*0.4; c=0.0003; 
        end
        pin=[p1,p2];
        E=p1/R;
        F=cross(cross(p1,p2),p1);
        F=F/norm(F);
  
        mu2=R*(E(1)*cos(x(:,3,i))+F(1)*sin(x(:,3,i)));
        mu1=-R*(E(2)*cos(x(:,3,i))+F(2)*sin(x(:,3,i)));
        nu=R*(E(3)*cos(x(:,3,i))+F(3)*sin(x(:,3,i)));

    ph = atan2(-mu1,mu2);
    th = acos(nu/R);
    
    path_lat=real((pi/2-th)*180/pi);
    path_long=real(ph*180/pi);
    geoshow(path_lat,path_long,'Color','w','LineWidth',2)
  
   axis([-.47,0.27,-1.57,-0.35])
end
set(gcf, 'Color', 'w')

