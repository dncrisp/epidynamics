function [isix,isiy,ampx,ampy,pks,locs] = isi_amp_onsetAnalysis(y,Fs,start,stop,thresh,diffthresh,ampwin)

% WARNING: this script does not filter your data. For best results, please
% filter your data before using this script.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%y: data vector
%
%Fs: Sampling Frequency. Used to recreate time vector
%
%start: Beginning time of the onsetseconds, with reference to the start of 
%       the data. Assumes first element of data vector occurs at time 1/FS.
%
%stop: The time at which you want to stop analyzing ISI and amplitude
%
%thresh: The threshold at which you want to detect spikes. For example, if
%        the maximum baseline amplitude is 0.1 units and the amplitude of
%        the spikes range between 0.4 and 0.6 units, then the threshold
%        should obviously be between 0.1 and 0.4 units to detect all the
%        spikes and none of the baseline noise.
%
%diffthresh: The temporal threshold at which you want to differentiate
%            spikes. For example, if the end of a seizure has clonic
%            bursting that lasts 0.6 seconds, and the clonic burst
%            complexes are separated by 1.4 seconds, then diffthresh should
%            be between 0.6 and 1.4 seconds. This would then treat the
%            clonic bursts as singular spikes, which is recommended. 
%
%ampwin: The window of time (in seconds) centered on a spike that searches
%        for amplitude. For example, assume a spike is detected at 100 
%        seconds. If ampwin = 4 seconds, then the script will look for the
%        smallest waveform value (gamma) in the range of 98 and 102. The 
%        peak-to-peak amplitude will then be calculated as: 
%        spike amplitude - gamma. A recommended value is 0.4 (or 400 ms). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS
%
%isiy: The inter-spike interval (ISI)
%
%isix: The times of the ISI
%
%ampy: The peak-to-peak amplitude of the spikes.
%
%ampx: The times of the spikes
%
%pks: spike peak absolute amplitude
%
%locs: spike peak times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning('off')

% Onset Analysis
time = (1:length(y))/Fs;
chosen_t = time(time > start & time < stop);
chosen_y = y(time > start & time < stop);
[pks, locs] = findpeaks(chosen_y, Fs, 'MinPeakHeight', thresh, 'MinPeakDistance',diffthresh);
locs2 = locs - locs(1);

% Time, ISI, and Amplitude
isiy = diff(locs2);
isix = locs2(1:end-1) + 0.1; % Assumption that last ISI point occurs 0.1 second away from seizure end.
ampx = locs2 + 0.1;
locs = locs + start;
for i = 1:length(pks)
    t_a = locs(i) - ampwin/2; 
    t_b = locs(i) + ampwin/2;
    temp_y = chosen_y(chosen_t > t_a & chosen_t < t_b);
    ampy(i) = pks(i) - min(temp_y);
end











