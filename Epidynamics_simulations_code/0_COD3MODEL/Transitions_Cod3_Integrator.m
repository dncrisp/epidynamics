% INTEGRATOR of Tranitions_Cod3_Model.m (Saggio et al. 2017, JMN)
% Transtions between classes implemented
% Marisa Saggio 2016

clear all;

%------- SETTINGS - INTEGRATION
options = odeset('MaxStep', 0.1);   % maximal time integration step

tstep=0.5; tmax=30000; tspan = [0:tstep:tmax];            % time span to integrate over
Nt=length(tspan);

x0in=[0;0;0;0;0];                           % initial conditions (must be a column)

%------- SETTINGS - MODEL
b = 1.0;                              % focus
R = 0.4;                            % radius in unfolding
c = 0.005;                           % velocity of slow variable
dstar =0.3 ;                          % threshold for slow variable inversion
N=1;
%------- PATH

% Path for BigTransition %axis([0,1700,-1.4,1.4])

c=0.002; cA=0.0001; cB=0.00012; tstep=0.5; tmax=18250; tspan = [0:tstep:tmax];

Ain=[0.2731,-0.05494,0.287];
Bin=[0.243,0.0461,0.3144];
Aend=[0.07337,-0.06485,-0.3878];
Bend=[-0.02792,-0.03676,-0.3973];


% c3 to c2
% tmax=5500; tspan = [0:tstep:tmax]; 
%  Ain=[0.2731,-0.05494,0.287];
%  Bin=[0.3331,0.074,0.2087];
%  Aend=[0.3524,0.05646,0.1806];
%  Bend=[0.3496,0.07955,0.1774];
%  c=0.0006;
%  cA=0.00008;
%  cB=0.000008;

% transition c10-c2
% Ain=[0.3131,-0.06743,0.2396];Bin=[0.29,0.06011,0.2689];
% Aend=[0.3448,0.02285,0.2014];Bend=[0.3496,0.07955,0.1774]; 
% c=0.0003;
% c2=0.00003;

% c10 to c3
% Ain=[0.3454,0.02484,0.2003];Bin=[0.29,0.06011,0.2689];
% Aend=[0.2731,-0.05494,0.287];Bend=[0.3331,0.074,0.2087];
% c=0.0001;
% cA=0.00001;cB=0.00001;
% tmax=35000; tspan = [0:tstep:tmax];
% axis_range=[0,30000,-1,0.8];

% c2 c2b
% Ain=[0.2731,-0.05494,0.287];
% %Bin=[0.2848,0.05851,0.2747];
% Bin=[0.243,0.0461,0.3144];
% Aend=[0.07337,-0.06485,-0.3878];
% Bend=[-0.02792,-0.03676,-0.3973];
% c=0.001
% c2=0.00008

%------- INTEGRATION
G=Ain/R;
H=cross(cross(Ain,Aend),Ain);
H=H/norm(H);

L=Bin/R;
M=cross(cross(Bin,Bend),Bin);
M=M/norm(M);

tic
[t,x]=ode45('Transitions_Cod3_Model',tspan, x0in, options, b, c,R,G,H,L,M,dstar,cA,cB,N); % evaluation of stochastic ODE
toc
%% ------- PLOTS

% plot curve followed in parameter space
openfig('FOCUS_r0p4/Unfolding_r0p4.fig')
hold on

    Au=R*(G'*cos(x(:,4))'+H'*sin(x(:,4))'); % A(u) and B(u)
    Bu=R*(L'*cos(x(:,5))'+M'*sin(x(:,5))');
    
    Eu=Au/R;                        % vectors for path from A(u) to B(u)
    Fu=cross(cross(Au,Bu),Au);
    for i=1:size(x(:,4),1)
    Fu(:,i)=Fu(:,i)/norm(Fu(:,i));
    end
   
    mu2=R*(Eu(1,:).*cos(x(:,3))'+Fu(1,:).*sin(x(:,3))'); 
    mu1=-R*(Eu(2,:).*cos(x(:,3))'+Fu(2,:).*sin(x(:,3))');
    nu=R*(Eu(3,:).*cos(x(:,3))'+Fu(3,:).*sin(x(:,3))');

plot3(mu2,-mu1,nu,'r','LineWidth',1.5)
set(gcf, 'Color', 'w')
hold off

%%

FigHandle = figure('Position', [100, 100, 70*(18/3), 70*3]);
newblue=[0, 113/255, 188/255];
newred=[216/255, 82/255, 24/255];
plot(tspan, x(:,1),'Color',newblue);
hold on
plot(tspan, x(:,3),'Color',newred);
plot(tspan, x(:,4));
plot(tspan, x(:,5));




