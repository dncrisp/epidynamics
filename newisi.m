function newisi(data,Fs)
% new isi det

% LOAD data first

isadaptive = 1;

y = [];         %#ok

% Fs = 512;
y = data(:,2);

locs    = [];   %#ok
peaks   = [];   %#ok
widths  = [];   %#ok
proms   = [];   %#ok

mpp     = 300.0;    % min peak prominence
mpd     = 0.0;      % min peak distance
mph     = 150.0;    % min peak height
mpw     = 0.0;      % min peak width
mxpw    = 0.1;      % max peak width
mth     = 0.0;      % min threshold 

sig_mpp = 0.3;
sig_mph = 0.30;

porder  = 7;
flength  = 41;   

windx  = [25 35];
yrange = 3.0;

yunfilt = y;
y = sgolayfilt(y,porder,flength);
    
if isadaptive == 1
    
    mpp     = 0.0;      % min peak prominence
    mpd     = 0.0;      % min peak distance
    mph     = 0.0;      % min peak height
    mpw     = 0.0;      % min peak width
    mxpw    = 0.2;     % max peak width
    mth     = 0.0;      % min threshold 
    
    mpp_d   = mpp;
    mpd_d   = mpd;
    mph_d   = mph;
    mpw_d   = mpw;
    mxpw_d  = mxpw;
    mth_d   = mth;

    mpp_u   = mpp;
    mpd_u   = mpd;
    mph_u   = mph;
    mpw_u   = mpw;
    mxpw_u  = mxpw;
    mth_u   = mth;
    
    
    [dpk,dlc,dwd,dpr] =  getpeaks(-y,Fs,mpp_d,mpd_d,mph_d,mpw_d,mxpw_d,mth_d);
    [upk,ulc,uwd,upr] =  getpeaks(y,Fs,mpp_u,mpd_u,mph_u,mpw_u,mxpw_u,mth_u);
    
    d_est_mph = positive(median(dpk) - sig_mph*std(dpk));
    d_est_mpp = positive(median(dpr) - sig_mpp*std(dpr));
    u_est_mph = positive(median(upk) - sig_mph*std(upk));
    u_est_mpp = positive(median(upr) - sig_mpp*std(upr));
    
    mpp_d = d_est_mpp;
    mph_d = d_est_mph;
    
    mpp_u = u_est_mpp;
    mph_u = u_est_mph;
    
    [dpk,dlc,dwd,dpr] =  getpeaks(-y,Fs,mpp_d,mpd_d,mph_d,mpw_d,mxpw_d,mth_d);
    [upk,ulc,uwd,upr] =  getpeaks( y,Fs,mpp_u,mpd_u,mph_u,mpw_u,mxpw_u,mth_u);
else
    
    [upk,ulc,uwd,upr] =  getpeaks( y,Fs,mpp,mpd,mph,mpw,mxpw,mth);
    [dpk,dlc,dwd,dpr] =  getpeaks(-y,Fs,mpp,mpd,mph,mpw,mxpw,mth);
end


% Plotting functionality -------------------------------

f = figure('units','normalized','position',[.1 .1 .7 .7]);

% Unfiltered data plot with overlaid filtered peaks
subplot(4,2,[1 2]);
plot(data(:,1),yunfilt,'k');
hold on
p2 = scatter(ulc,upk,'o');
p3 = scatter(dlc,-dpk,'o');
legend([p2,p3],'Up Spikes','Down Spikes');
hold off
 xlim(windx);
 yiqr = iqr(data(:,2));
 ylim([-yiqr*yrange yiqr*yrange]);
xlabel('Unfiltered Data with Filtered Peaks');

tstring = sprintf('Peak finder parameter tuning (down): mpp=%.1f, mpd=%.1f, mph=%.1f, mpw=%.1f, w/sgfilt',...
    mpp_d, mpd_d, mph_d, mpw_d);
% tstring = sprintf('Peak finder parameter tuning (down): mpp=%.1f, mpd=%.1f, mph=%.1f, mpw=%.1f, w/sgfilt',...
%     mpp, mpd, mph, mpw);
title(tstring);

% Filtered data plot with 
subplot(4,2,[3 4]);
plot(data(:,1),y,'k');
hold on
scatter(ulc,upk,'o');
scatter(dlc,-dpk,'o');
 xlim(windx);
 yiqr = iqr(data(:,2));
 ylim([-yiqr*yrange yiqr*yrange]);
hold off
xlabel(sprintf('Filtered: sgolay(%d,%d)',porder,flength));


subplot(4,2,5);
histogram(upk);
hold on
histogram(-dpk);
title('Histogram of Peak Amplitude');
hold off

piu = diff(ulc);
pid = diff(dlc);

subplot(4,2,6);
histogram(piu);
hold on
histogram(pid);
title('Histogram of Peak Interval');
hold off
xlim([0 1]);

subplot(4,2,7);
histogram(uwd);
hold on
title('Histogram of Peak Widths');
histogram(dwd);
hold off
xlim([0 0.5]);

subplot(4,2,8);
histogram(upr);
hold on
histogram(dpr);
title('Histogram of Peak Prominence');
hold off

szend = dlc(end);
szneg = dlc(:,1) - szend;
szrevx = abs(flipud(szneg));
szrevy = szrevx;

isix = szrevx(1:end-1);
isiy = diff(szrevy);

figure;



scatter(isix,isiy);

end



function [locs, peaks, widths, proms] = getpeaks(y,Fs,mpp,mpd,mph,mpw,mxpw,mth)

[locs,peaks,widths,proms] = findpeaks(y,Fs,...
    'MinPeakProminence',mpp,...
    'MinPeakDistance',mpd,...
    'MinPeakHeight',mph,...
    'MinPeakWidth',mpw,...
    'MaxPeakWidth',mxpw,...
    'Threshold',mth,...
    'WidthReference','halfheight',...
    'Annotate','extents');

end


function value = positive(x)

    if x < 0
        value = 0;
    else 
        value = x;
    end

end    
    
% yy = smooth(isix,isiy,0.005,'rloess');
% %yy = sgolayfilt(isiy,7,101);
% 
% figure;
% plot(isix,yy);
