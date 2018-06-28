%% Simulations for Fig_S11
% Marisa Saggio
%
% Class SN/supH, with three different paths giving different trends for the
% frequency behavior
% The model used is described in Saggio et al. 2017, JMN

%%
clear all;

%------- SETTINGS - INTEGRATION
options = odeset('MaxStep', 0.1);   % maximal time integration step
dt=0.1;
tspan=[0:0.5:20000];              % time span to integrate over
x=zeros(length(tspan),3,4);        % initial conditions 

%------- SETTINGS MODEL
b = 1.0;                            % deg. TB singularity, focus case
R = 0.4;                            % radius in unfolding
c = 0.001;                          % velocity of slow variable
dstar =0.3 ;                        % threshold for slow variable inversion
N=1;                                % if N=1, upper branch of equilibria is resting state

%% INTEGRATION

for i=1:3
    
    % sets path parameters. p1 is called A in the paper, p2 is called B.
    switch i
        case 1 % constant frequency
            p1=[0.7753,-0.05421,0.6293]*0.4;p2=[0.8564,0.1929,0.4788]*0.4;c = 0.00025;
            class=sprintf('c2s - SN/SH');
            range_axis=[4000,6700,-1.2,0.7];
        case 2 % increasing frequency
            p1=[0.8041,0.1563,0.5736]*0.4;p2=[0.8697,0.1974,0.4524]*0.4; c=0.001;
            class=sprintf('c3s - increasing frequency');
            range_axis=[3700,5200,-1.2,0.7];
        case 3 % decreasing frequency
            p2=[0.3249,0.07129,0.2221];tspan=[0:0.5:20000];
            p1=[0.3102,-0.04369,0.2488]; c=0.0002;
            class=sprintf('c3s - decreasing frequency');
            range_axis=[4900,7500,-1,0.8];
    end
    
    % computes path vectors E,F
    pin=[p1,p2];
    E=p1/R;
    F=cross(cross(p1,p2),p1);
    F=F/norm(F);
    
    % integrates equations
    [t,x(:,:,i)]=ode45('Cod3_Model',tspan, x0, options, b, c, R, dstar,E,F,N);
end

%% PLOT TIMESERIES
figure
for i=1:3
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

% plot path on Amplitude map

openfig('../../0_COD3MODEL/Figures/Ampl_Equazim.fig')
hold on
topolegend=[1,90,0];
axesm('MapProjection','eqaazim','Grid','on')
for i=1:3
    switch i
        case 1
            p1=[0.7753,-0.05421,0.6293]*0.4;p2=[0.8564,0.1929,0.4788]*0.4;
        case 2
            p1=[0.8041,0.1563,0.5736]*0.4;p2=[0.8697,0.1974,0.4524]*0.4;
        case 3
            p2=[0.3249,0.07129,0.2221];p1=[0.3102,-0.04369,0.2488];
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
    % plot
    geoshow(path_lat,path_long,'Color','w','LineWidth',2)
end
set(gcf, 'Color', 'w')

% plot path on Frequency map

openfig('../../0_COD3MODEL/Figures/Freq_Equazim.fig')
hold on
topolegend=[1,90,0];
axesm('MapProjection','eqaazim','Grid','on')
for i=1:3
    switch i
        case 1
            p1=[0.7753,-0.05421,0.6293]*0.4;p2=[0.8564,0.1929,0.4788]*0.4;
        case 2
            p1=[0.8041,0.1563,0.5736]*0.4;p2=[0.8697,0.1974,0.4524]*0.4;
        case 3
            p2=[0.3249,0.07129,0.2221];p1=[0.3102,-0.04369,0.2488];
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
    % plot
    geoshow(path_lat,path_long,'Color','w','LineWidth',2)
end
set(gcf, 'Color', 'w')

