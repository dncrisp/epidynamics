% plot for Fig S6, B-C
figure

subplot(2,1,1)
load('../Status_beforeSH.mat') % this is the 4th simulation of RepRep2=4
tmin_plot=3000/dt;
tmax_plot=6000/dt;
plot(tspan(tmin_plot:tmax_plot),X(tmin_plot:tmax_plot,1))

subplot(2,1,2)
load('Status_belowSH.mat') % this is the 10th simulation of RepRep2=4
tmin_plot=5000/dt;
tmax_plot=10000/dt;
plot(tspan(tmin_plot:tmax_plot),X(tmin_plot:tmax_plot,1))