% Produces simulations from the model in Saggio et al. 2017 JMN
%
% Uses the functions:   Parametrization_VectorsGreatCircle.m
%                       Cod3_Model.m
%                           Parametrization_GreatCircle.m
%                           eval_resting_state_cartesian.m
%
% Plotting the results uses the figure FOCUS_r0p4/Unfolding_r0p4.fig
%
% Examples of settings for different classes are given (commented) in 'Path
% Settings'

clear all;

%% General Settings
%------- SETTINGS - INTEGRATION
options = odeset('MaxStep', 0.1);   % maximal time integration step
tspan = [0:0.5: 3000 ];               % time span to integrate over
x0=[0;0;0];                         % initial conditions (must be a column)

%------- SETTINGS - MODEL
b = 1.0;                            % focus
R = 0.4;                            % radius in unfolding
c = 0.0003;                           % velocity of slow variable
dstar =0.3;                        % threshold for slow variable inversion
N=1;                                % with N=1 the upper branch of equilibria is the resting state

%% Path Settings
% initial (on bif curve which stops the oscillations) point A and final
% point B (on bif curve which starts the oscillations) in cartesian 
% parameter space for the system trajectory (x,y,z) in figure of unfolding (meaning (mu2,-mu1,nu))

% Setting A,B
% SMALL LC region

%A=[0.2649,-0.05246,0.2951];B=[0.2688,0.05363,0.2914];%c = 0.002;tspan=[0:0.1:100 ]; % SN/SN
A=[0.3448,0.02285,0.2014];B=[0.3351,0.07465,0.2053]; c=0.0001;tspan = [0:0.5:10000 ];   % SN/SH  
%A=[0.2552,-0.0637,0.3014];B=[0.3496,0.0795,0.1774];c=0.0004;tspan =[0:0.5:10000 ]; % SN/SupH
%A=[0.3448,0.0228,0.2014]; B=[0.3118,0.0670,0.2415];c = 0.00005; tspan =[0:0.05:20000 ]; % SupH/SH
% A=[0.3131,-0.06743,0.2396];B=[0.3163,0.06846,0.2351]; %c=0.00004;  tspan = [0:0.05:50000 ]; % SupH/SupH
% A=[0.3085,-0.06595,0.2459];B=[0.3207,0.06989,0.2287]; %c=0.000015; tspan =[0:0.5:150000]; % SupH/SupH bis

% BIG LC region

%A=[0.3216,0.0454,-0.2335];B=[0.285,0.05855,-0.2745]; c=0.004; tspan =[0:0.05:400 ];% SN/SHb
%c=0.004; dstar=0.3; B=[0.3883;0.03687;-0.2521];B=[0.184;0.02903;-0.354];tspan = [0:0.05:1400 ];% SubH/SH
%A=[0.3216,0.0454,-0.2335];B=[0.106,0.005238,-0.3857]; % SubH/SH bis
%A=[0.1871,-0.02512,-0.3526];B=[0.2081,-0.01412,-0.3413];c=0.008;tspan = [0:0.05:2000 ]; dstar=0.1% class SN/FLC 
%B=[-0.01301,-0.03242,-0.3985];A=[0.04098,-0.07373,-0.391];    % subH/FLC

%% Compute path and integrate

[E,F]=Parametrization_VectorsGreatCircle(A,B,R); % computes paths

[t,x]=ode45('Cod3_Model',tspan, x0, options, b, c, R, dstar,E,F,N); % evaluation of ODE

%% ------- PLOTS

% plot path followed in parameter space in red
openfig('FOCUS_r0p4/Unfolding_r0p4.fig')
hold on
[mu2, mu1, nu]=Parametrization_GreatCircle(E,F,R,x(:,3));
plot3(mu2,-mu1,nu,'r','LineWidth',1.5)
hold off

% plot timeseries
FigHandle = figure('Position', [100, 100, 850, 200]);
str=sprintf('A = %.4f,%.4f,%.4f;  B = %.4f,%.4f,%.4f, c = %.5f, dstar = %.2f', A(1), A(2), A(3), B(1), B(2), B(3), c, dstar);
title(str)
subplot(1,3,2:3)
newblue=[0, 113/255, 188/255];
newred=[216/255, 82/255, 24/255];
plot(tspan, x(:,1),'Color',newblue);
hold on
plot(tspan, x(:,3),'Color',newred);
hold off
title(str)
xlabel('time [s]')
legend('x','z')

% plot state space
subplot(1,3,1)
plot3(x(:,1),x(:,2),x(:,3),'Color',newblue)
xlabel('x')
ylabel('y')
zlabel('z')
set(gcf, 'Color', 'w')
  
