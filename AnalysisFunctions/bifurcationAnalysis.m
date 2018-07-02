% Data analysis code for the paper Epidynamics: Navigating the map of
% seizure dynamics
% This code has been set up to analyze the data for Figure S1A

clear all; close all; clc;
load('...\Fig_S1\Fig_S1A.mat') % Load Data
addpath('...\AnalysisFunctions') % Include path to functions
Fs = 1/mean(diff(data(2,:)));

% Filter Settings
Fstop1 = 1;
Fpass1 = 2;
Fpass2 = 55;
Fstop2 = 70;
Astop1 = 98;
Apass  = 0.5;
Astop2 = 100;
d = designfilt('bandpassiir', ...
  'StopbandFrequency1',Fstop1,'PassbandFrequency1', Fpass1, ...
  'PassbandFrequency2',Fpass2,'StopbandFrequency2', Fstop2, ...
  'StopbandAttenuation1',Astop1,'PassbandRipple', Apass, ...
  'StopbandAttenuation2',Astop2, ...
  'DesignMethod','butter','SampleRate', Fs);

% Display Raw Data
figure(1)
t = data(2,:);
plot(t,data(1,:), 'black')
title('Raw Data')
xlabel('Time (s)')

% Display Filtered Data
figure(2)
temp = filtfilt(d, data(1,:));
%temp = sgolayfilt(temp, 7, 41); %Optional filter
r = 20;
dataf = decimate(temp, r); 
t2 = (0:(length(dataf)-1))/(Fs/r);
plot(t2,dataf, 'black')
ylim([-3000 5000])
title('Filtered Data')
xlabel('Time (s)')

% Select Onset Window
a = 139; % Start of onset window to be analyzed
b = 145; % End of onset window to be analyzed
tonset = t2(t2 > a & t2 < b);
donset = dataf(t2 > a & t2 < b);

% Initialize fit equations
ft_log = fittype('a*log(x) + b', 'coeff', {'a', 'b'});
ft_con = fittype('a*x*0 + b', 'coeff', {'a', 'b'});
ft_sqr = fittype('a*sqrt(x) + b', 'coeff', {'a', 'b'});

% Spike Extraction
ampThresh = 50; % MinPeakHeight
timeThresh = 0.009; %MinPeakDistance
ampWin = 0.4; % Window for amplitude analysis

% For Onset Analysis:
[isix,isiy,ampx,ampy,pks,locs] = isi_amp_onsetAnalysis(dataf,Fs/r,a,b,ampThresh,timeThresh,ampWin);
% For Offset Analysis
%[isix,isiy,ampx,ampy,pks,locs] = isi_amp_offsetAnalysis(dataf,Fs/r,a,b,ampThresh,timeThresh,ampWin);
isix = isix';
isiy = isiy';
ampx = ampx';
ampy = ampy';

% Plot Onset Window
figure(3)
subplot(3,1,1)
plot(tonset, donset, 'black')
xlim([ a-0.5 b+0.5])
title('ISI & Amp Analysis')

% Select seizure spikes
% This should be used to eliminate false spikes
goodSpikes = [11,12,13,15:length(locs)];
locs = locs(goodSpikes);
pks = pks(goodSpikes);

ampx = ampx(goodSpikes);
ampx = ampx - min(ampx)+ 0.1;
ampy = ampy(goodSpikes);

% Used to index reasonable ISI data points
goodISI = [11,12,15:length(locs)];
isix = isix(goodISI);
isix = isix - min(isix) + 0.1;
isiy = isiy(goodISI);


% Plot used spikes on filtered data
hold on
plot(locs,pks,'o', 'Color', 'red')
hold off

% Plot ISI and best fit equation
subplot(3,1,2)
plot(isix, isiy, 'o', 'Color', 'red')
hold on
[eqfit1,gof1,out1] = fit(isix,isiy,ft_con);
x=linspace(isix(1),isix(end),100);
plot(x, eqfit1.a*log(x)*0 + eqfit1.b, 'black');
xlim([min(isix)-0.5 max(isix)+0.5])
hold off

% Plot Amplitude and best fit equation
subplot(3,1,3)
plot(ampx, ampy, 'o', 'Color', 'red')
hold on
[eqfit2,gof2,out2] = fit(ampx,ampy,ft_con);
plot(ampx, eqfit2.a*sqrt(ampx)*0 + eqfit2.b, 'black');
xlim([min(ampx)-0.5 max(ampx)+0.5])
hold off

% GOF testing
clc
[~,gof_isi_sqr,~] = fit(isix,isiy,ft_sqr)
[~,gof_isi_con,~] = fit(isix,isiy,ft_con)
[~,gof_amp_sqr,~] = fit(ampx,ampy,ft_sqr)
[~,gof_amp_con,~] = fit(ampx,ampy,ft_con)



