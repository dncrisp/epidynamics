%% Fig 5, H-I

%% Status epilepticus

load('Status_def.mat') % this is the first simulation of RepRep2=1

fs_exp = 512;   %[fs]=1/s %system sampling frequency
dt_exp = 1 / fs_exp;  %[dt]=s
dt=dt_exp; tmax= 20000;  tspan = (0 : dt : tmax)';    % time span to integrate over
xinit=[0;0;0;0;0];            % initial conditions (must be a column)

b = 1.0; R = 0.4; dstar =0.3 ; N_equil=1;  

%noise
sigma=0.0005*[10,10,1,5,5];
dim = [1,2,3,4,5]; %additive correlated noise to variable in dim

%path
c = 0.0001;cA = 0;cB=0; 
Ain=[0.3483,0.03698,0.1931];
Bin=[0.3331,0.074,0.2087];
Aend=[0.279,0.2187,0.1854];
Bend=Aend;

% vectors G,H,M,L to describe path
G=Ain/R; H=cross(cross(Ain,Aend),Ain); H=H/norm(H);
L=Bin/R; M=cross(cross(Bin,Bend),Bin); M=M/norm(M);

% plot curve followed in parameter space
openfig('../0_COD3MODEL/FOCUS_r0p4/Unfolding_r0p4.fig')
hold on

    Au=R*(G'*cos(X(:,4))'+H'*sin(X(:,4))'); % A(u) and B(u)
    Bu=R*(L'*cos(X(:,5))'+M'*sin(X(:,5))');
    
    Eu=Au/R;                        % vectors for path from A(u) to B(u)
    Fu=cross(cross(Au,Bu),Au);
    for i=1:size(X(:,4),1)
    Fu(:,i)=Fu(:,i)/norm(Fu(:,i));
    end
   
    mu2=R*(Eu(1,:).*cos(X(:,3))'+Fu(1,:).*sin(X(:,3))'); 
    mu1=-R*(Eu(2,:).*cos(X(:,3))'+Fu(2,:).*sin(X(:,3))');
    nu=R*(Eu(3,:).*cos(X(:,3))'+Fu(3,:).*sin(X(:,3))');

    tmin_plot=3510/dt;
tmax_plot=6503/dt;

plot3(mu2(tmin_plot:tmax_plot),-mu1(tmin_plot:tmax_plot),nu(tmin_plot:tmax_plot),'k','LineWidth',1.5)
%plot3(mu2,-mu1,nu,'r','LineWidth',1.5)
set(gcf, 'Color', 'w')
hold off
%% plot flat

ph = atan2(-mu1,mu2);
th = acos(nu./(sqrt(mu2.^2+mu1.^2+nu.^2)));

% plot freq

openfig('../0_COD3MODEL/Figures/Freq_Equazim.fig')
hold on
topolegend=[1,90,0];
axesm('MapProjection','eqaazim','Grid','on')
   
    path_lat=(pi/2-th)*180/pi;
    path_long=ph*180/pi;
    geoshow(path_lat(tmin_plot:tmax_plot),path_long(tmin_plot:tmax_plot),'Color','r','LineWidth',3)
    geoshow(path_lat(tmin_plot:tmax_plot),path_long(tmin_plot:tmax_plot),'Color','w','LineWidth',2)
    plotm(path_lat(tmin_plot),path_long(tmin_plot),'^r','LineWidth',0.5,'MarkerFaceColor','r','MarkerSize',6)
    plotm(path_lat(tmax_plot),path_long(tmax_plot),'xr','LineWidth',2,'MarkerSize',8)
   
    title('Frequency')
    xlabel('\phi')
    ylabel('\theta')
    
    % plot ampl

openfig('../0_COD3MODEL/Figures/Ampl_Equazim.fig')
hold on
topolegend=[1,90,0];
axesm('MapProjection','eqaazim','Grid','on')
   
    path_lat=(pi/2-th)*180/pi;
    path_long=ph*180/pi;
    geoshow(path_lat(tmin_plot:tmax_plot),path_long(tmin_plot:tmax_plot),'Color','w','LineWidth',2)
    plotm(path_lat(tmin_plot),path_long(tmin_plot),'^r','LineWidth',0.5,'MarkerFaceColor','r','MarkerSize',6)
    plotm(path_lat(tmax_plot),path_long(tmax_plot),'xr','LineWidth',2,'MarkerSize',8)

    title('Amplitude')
    xlabel('\phi')
    ylabel('\theta')
    
%% c3-c2

load('transition_within_seizure.mat')

fs_exp = 512;                 %[fs]=1/s %system sampling frequency
dt_exp = 1 / fs_exp;               %[dt]=s
dt=dt_exp; tmax= 10000; tspan = (0 : dt : tmax)';    % time span to integrate over
xinit=[0;0;0;0;0];            % initial conditions (must be a column)

%------- SETTINGS - MODEL
b = 1.0; R = 0.4; dstar =0.3 ; N_equil=1;                          % resting state is upper branch of equilibria

%noise
sigma=0.0005*[10,10,1,4,4];
dim = [1,2,3,4,5]; %additive correlated noise to variable in dim

% % settings used to have transition in between def
c = 0.0002; cA = 0.00005; cB=0.000005; 
 Ain=[0.2731,-0.05494,0.287];
 Bin=[0.3331,0.074,0.2087];
 Aend=[0.3524,0.05646,0.1806];
 Bend=[0.3496,0.07955,0.1774];
 
 % vectors G,H,M,L to describe path
G=Ain/R; H=cross(cross(Ain,Aend),Ain); H=H/norm(H);
L=Bin/R; M=cross(cross(Bin,Bend),Bin); M=M/norm(M);

% plot curve followed in parameter space
openfig('../0_COD3MODEL/FOCUS_r0p4/Unfolding_r0p4.fig')
hold on

    Au=R*(G'*cos(X(:,4))'+H'*sin(X(:,4))'); % A(u) and B(u)
    Bu=R*(L'*cos(X(:,5))'+M'*sin(X(:,5))');
    
    Eu=Au/R;                        % vectors for path from A(u) to B(u)
    Fu=cross(cross(Au,Bu),Au);
    for i=1:size(X(:,4),1)
    Fu(:,i)=Fu(:,i)/norm(Fu(:,i));
    end
   
    mu2=R*(Eu(1,:).*cos(X(:,3))'+Fu(1,:).*sin(X(:,3))'); 
    mu1=-R*(Eu(2,:).*cos(X(:,3))'+Fu(2,:).*sin(X(:,3))');
    nu=R*(Eu(3,:).*cos(X(:,3))'+Fu(3,:).*sin(X(:,3))');

    tmin_plot=1450000; tmax_plot=2750000;
plot3(mu2(tmin_plot:tmax_plot),-mu1(tmin_plot:tmax_plot),nu(tmin_plot:tmax_plot),'r','LineWidth',1.5)
%plot3(mu2,-mu1,nu,'r','LineWidth',1.5)
set(gcf, 'Color', 'w')

%% plot flat

ph = atan2(-mu1,mu2);
th = acos(nu./(sqrt(mu2.^2+mu1.^2+nu.^2)));

% plot freq

openfig('../0_COD3MODEL/Figures/Freq_Equazim.fig')
hold on
topolegend=[1,90,0];
axesm('MapProjection','eqaazim','Grid','on')
   
    path_lat=(pi/2-th)*180/pi;
    path_long=ph*180/pi;
    geoshow(path_lat(tmin_plot:tmax_plot),path_long(tmin_plot:tmax_plot),'Color','w','LineWidth',2)
    plotm(path_lat(tmin_plot),path_long(tmin_plot),'^r','LineWidth',0.5,'MarkerFaceColor','r','MarkerSize',6)
    plotm(path_lat(tmax_plot),path_long(tmax_plot),'xr','LineWidth',2,'MarkerSize',8)

    title('Frequency')
    xlabel('\phi')
    ylabel('\theta')
    
    % plot ampl

openfig('../0_COD3MODEL/Figures/Ampl_Equazim.fig')
hold on
topolegend=[1,90,0];
axesm('MapProjection','eqaazim','Grid','on')
   
    path_lat=(pi/2-th)*180/pi;
    path_long=ph*180/pi;
    geoshow(path_lat(tmin_plot:tmax_plot),path_long(tmin_plot:tmax_plot),'Color','w','LineWidth',2)
    plotm(path_lat(tmin_plot),path_long(tmin_plot),'^r','LineWidth',0.5,'MarkerFaceColor','r','MarkerSize',6)
    plotm(path_lat(tmax_plot),path_long(tmax_plot),'xr','LineWidth',2,'MarkerSize',8)

    title('Amplitude')
    xlabel('\phi')
    ylabel('\theta')
    