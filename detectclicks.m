function [click,fs,newthreshold] = detectclicks (file,threshold,SETTINGS)
% [click,fs,threshold] = detectclicks(file,threshold,SETTINGS)
% Finds the samples corresponding to the peak of clicks exceeding threshold
% 
% Input:
%   file:       character string containing path and file name (including
%               extension) of wave file compatible with wavexread
%   threshold:  Detection threshold (relative to 1). Can be left blank to
%               implement an adaptive threshold
% SETTINGS is a structure with at least the following fields:
%  CH:          Channel to search for clicks
%  detect_filt: High-pass [HP] or bandpass [HP LP] filter (in Hz)
%  detect_bl:   Click detector blanking time (in seconds)
%  clickdur:    Expected click duration (in seconds)
%
% [click,fs,newthreshold] = detectclicks(file,[],SETTINGS)
% with no threshold given, click detector will use an adaptive threshold
% based on noise of individual chunks
%
% Output:
%   click:      List of samples corresponding to the peak of the envelope
%               for each click exceeding threshold
%   fs:         Sample rate of wave file
%   threshold:  Will return either the constant threshold used, or, in the
%               case of no given threshold, a list of all adaptive
%               threshold values
%
% F. H. Jensen, 2013 (frants.jensen@gmail.com)

% Load settings
if nargin<3,
    CH      = 3;      % Default channel to perform click detection on
    HP      = 20e3 ;  % Default high-pass filter
    blank   = 0.002 ; % Default blanking time
    clickdur = 1e-5 ; % Default expected click duration   
else
    CH      = SETTINGS.CH ;
    HP      = SETTINGS.detect_filt ;
    blank   = SETTINGS.detect_bl ;
    clickdur = SETTINGS.clickdur ;
end

% Check for threshold
if ~exist('threshold')
    nothr = 1;                                  % Activate adaptive threshold  
    threshold = [];                             % Adaptive threshold vector
elseif isempty(threshold) | threshold ==0,
    nothr = 1;                                  % Activate adaptive threshold
    threshold = [];                             % Adaptive threshold vector
else
    nothr = 0;
end

% Prepare click detector

% Find length and sample rate of wav file
file_size = wavexread(file,'size');             % Find total size of sound
file_size = file_size(1);                       % Discard number of channels
[sig,fs] = wavexread ( file , [1 5] );          % Find fs

% Design filter
if length(HP)==2
    [B,A] = butter(4,HP/(fs/2)) ;               % Bandpass filter
else
    [B,A] = butter(4,HP/(fs/2),'high');         % Highpass filter
end

% Constants:
chunk = 200000 ;        % Number of samples to be read at once - 50k optimized
jump  = blank*fs ;      % Blanking time in samples
N1    = jump ;          % Define default first sample to be analysed
dur   = ceil(2*clickdur*fs); % Samples to search for maximum envelope detection
n     = 0;              % Detected click counter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following while-loop will run until the end of the file
hwaitbar = waitbar(0,' Detecting clicks ');
while N1 < file_size-jump

    % Define part of sound file to read
    if N1+chunk-1 > file_size, % If close to end of file, take last bit
        [sig,fs] = wavexread (file , [N1 file_size]);
    else % Otherwise read in full chunk
        [sig,fs] = wavexread (file , [N1 N1+chunk-1]);
    end
    
    % Clear memory of unnecessary data    
    sig = sig(:,CH); 

    % Warn if data are clipped
    if max(abs(sig))>=0.99
        disp ('warning - input data might be clipped' );
        disp( ['progress: N1 = ' num2str(N1) ] )
    end

    % Filter signal
    sig = sig-mean(sig);
    sig = filter(B,A,sig);
    env = abs(hilbert(sig));

    % Adaptive threshold technique - from Mark Johnson rainbow click detector (rainbow.m)
    if nothr
        thr = raylinv(0.9999,raylfit(env)) ;
        thr = thr*4 ; % No minimum threshold, no maximum threshold, threshold factor 4 (conservative)
        threshold = [threshold ; thr];
    else
        thr = threshold ;
    end

    % Find clicks exceeding threshold
    k = find(env>thr);
    if ~isempty(k)
        firstsample = k([1 ; 1+find(diff(k)>jump)]);
        for int = 1:length(firstsample)
            if k(int)+dur>length(env)
                continue
            end
            [temp,peak] = max(env(k(int)+[0:dur]));
            n = n+1;
            click(n) = (N1-1) + (firstsample(int)-1) + peak;
        end
    end
    
    waitbar(N1/file_size,hwaitbar)
    
    % Once end of chunk has been reached, proceed with next chunk:
    N1  = N1 + chunk - dur;
     
end

% Finish and close waitbar
waitbar(1,hwaitbar)
close(hwaitbar)

% Dismiss clicks that are too close to end of file
if ~isempty(click)
    click=click(find(click<file_size-jump));
end
    

% Evaluate click detection
disp(['Click Detection Complete for file ' file])
if ~exist('click'), error('No clicks detected in file - adjust threshold?'),end
disp([num2str(max([N1-jump click(end)])) ' samples out of ' num2str(file_size) ' analysed and ' num2str(n) ' clicks detected.']) 

newthreshold = threshold ;