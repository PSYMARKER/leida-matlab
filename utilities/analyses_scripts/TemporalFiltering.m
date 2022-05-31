function signal_filt = TemporalFiltering(signal,flp,fhi,tr)

% Temporal Filtering Settings
fnq = 1/(2*tr);                 % Nyquist frequency
Wn = [fhi/fnq flp/fnq];         % butterworth bandpass non-dimensional frequency
k = 5;                          % 5th order butterworth filter
[bfilt,afilt] = butter(k,Wn);   % construct the filter

% Number of brain areas
N = size(signal,1);
% Variable to store the filtered signal
signal_filt = zeros(N,size(signal,2));

for seed = 1:N
    signal_seed = signal(seed,:);
    % Apply the temporal filter to the signal of region=seed
    signal_filt(seed,:) = filtfilt(bfilt,afilt,signal_seed);
end