function [time,maxenv] = findthreshold(file,SETTINGS)
% [time,maxenv]=findthreshold(file,SETTINGS);
%
% Goes through file and extracts peak of hilbert envelope
% for each block of 10 ms, then returns a vector of time 
% and maximum envelope value
%
% SETTINGS is a structure with fields:
%  CH:              Channel to search for clicks
%  detect_filt:     High-pass or bandpass filter (in Hz)
%
% Requires multi-channel wav-file (compatible with wavexread) as input
%
% F. H. Jensen, 2013 (frants.jensen@gmail.com)

if nargin<2,
    CH = 3;     % Default channel to perform click detection on
    HP = 20e3 ; % Default high-pass filter
else
    CH = SETTINGS.CH ;
    HP = SETTINGS.detect_filt ;
end

% Find file size and sample rate
[file_size] = wavexread ( file , 'size' );
[x fs] = wavexread(file , [1 2]);

% Define chunk size to be loaded at a time
chunk = round(0.01*fs); % Default block size is 10 ms

% Design filter
if length(HP)==2
    [B,A] = butter(4,HP/(fs/2)) ;
else
    [B,A] = butter(4,HP/(fs/2),'high');
end

% Find maximum envelope point for each chunk
maxenv = zeros(floor(file_size(1)/chunk),1);
for i=1:floor(file_size(1)/chunk)
    N = [(i-1)*chunk+1 min([i*chunk file_size(1)])];
    [x fs] = wavexread ( file , N );
    xf = filter(B,A,x(:,CH)-mean(x(:,CH)));
    xh = abs(hilbert(xf));
    maxenv(i) = max(xh);
end

time = [1:length(maxenv)]/(fs/chunk);