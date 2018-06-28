%% Fig 5, A-F

%%
% load experimental data
load('EXP_DATA_switchingBifurcations.mat')

%% SN/supH - SN/SH transition within seizure

load('transition_within_seizure.mat')

figure
set(gcf, 'Color', 'w')

fr=HUP83_phaseII_fs;
points=length(HUP83_phaseII);
tmax=points/fr;
dt=1/fr;
tspan=0:dt:(tmax-dt);

subplot(3,1,1)
plot(tspan,HUP83_phaseII)
ylabel('AC recording')
xlabel('time [s]')

tmin_plot=1450000; tmax_plot=2750000;

points=length(X(:,1));
tmax=points/fr;
dt=1/fr;
tspan=0:dt:(tmax-dt);
subplot(3,1,2)
Y=highpassfilter_mine(X(:,1),0.03,1);
%Y=highpassfilter_mine(X(:,1),0.16,5); %tim filter
plot(tspan(tmin_plot:tmax_plot),Y(tmin_plot:tmax_plot))
ylabel('High-pass filtered simulated x')

subplot(3,1,3)
plot(tspan(tmin_plot:tmax_plot),X(tmin_plot:tmax_plot))
ylabel('Simulated x')

%% Status epilepticus

load('Status_def.mat') % this is the first simulation of RepRep2=1

figure
set(gcf, 'Color', 'w')

fr=Status_I002_P005_D01_fs;
points=length(Status_I002_P005_D01);
tmax=points/fr;
dt=1/fr;
tspan=0:dt:(tmax-dt);

tmin_plot=717/dt;
tmax_plot=1093/dt;
subplot(3,1,1)
plot(tspan(tmin_plot:tmax_plot),Status_I002_P005_D01(tmin_plot:tmax_plot))
ylabel('AC recording')
xlabel('time [s]')

%tmin_plot=1000000; tmax_plot=3370000;
tmin_plot=3501/dt;
tmax_plot=6501/dt;

points=length(X(:,1));
tmax=points/fr;
dt=1/fr;
tspan=0:dt:(tmax-dt);

subplot(3,1,2)
Y=highpassfilter_mine(X(:,1),0.03,1);
plot(tspan(tmin_plot:tmax_plot),Y(tmin_plot:tmax_plot))
ylabel('High-pass filtered simulated x')

subplot(3,1,3)
plot(tspan(tmin_plot:tmax_plot),X(tmin_plot:tmax_plot))
ylabel('Simulated x')