function S = toolbox_settings()
% Configure settings for click parameter toolbox
%
% Change this file before using click parameter toolbox
% - Adjusting these parameters correctly is critical for good performance
%   of click parameter toolbox
%
% F. H. Jensen, 2013 (frants.jensen@gmail.com)

% Setup file name
setupname = 'Snub';

% Define directory for output data files (click parameters) and prefix
datadir     = 'C:\Matlab_Scripts\biosonar_array_toolbox\Sousa\' ; % Click parameter files will be saved here, and click detections will be saved in subfolder
clickprefix = 'snubfin' ;

% Recording equipment
gain        = 40 ;          % dB amplification of system for this file
dbmaxvoltage= 0 ;          % Peak voltage of ADC in dB; e.g. 20log10(5V peak) = 14 dB re 1V
sens        = -218 ;        % Hydrophone nominal sensitivity, dB re 1V/muPa
sens_corr   = [0 0 0 0];  % Sensitivity correction for each hydrophone (near centroid frequency of species)
clip_level  = dbmaxvoltage - (sens+sens_corr) - gain  ; % Hydrophone peak clip level for each channel

% Array configuration
channels = [1 2 3 4] ;    % Channels used in click localization - must match length of depth vector
channel_center = 3;         % Channel used for click detection (center channel)
channel_depths = -0.95+[0:-0.90:-2.7]';  % Hydrophone depths (must be negative)

% Environment
soundspeed  = 1529; % Sound speed - use soundspeed.m to convert from temp, salinity and depth
temperature = 24 ;

% Click detector settings
detect_filt = [10e3] ;      % Click detector filter (high-pass [FH] or bandpass [FH FL])
detect_bl   = 0.002 ;       % Blanking time to skip after click detected (2 ms = 500 clicks per sec)
detect_clickdur = 10e-6 ;   % Expected click duration in seconds

% Click localization settings
localize_method = 'correlation' ; % Method of finding time difference: 
                                % Acceptable input {'peak','correlation','6db'}
                                % TOADs either peak envelope, cross-correlation (best for most things), 
                                % or first -6 dB shoulder (good for NBHF)

% Click parameter settings
NFFT            = 32 ;          % FFT size (number of samples extracted per click)
FFT_interp      = 10 ;         % FFT interpolation
centroid        = 90 ;          % Approximate on-axis centroid frequency for estimating absorption loss (dB/m)
filt            = [4000] ;      % Filter to apply during analysis, Hz (either HP or Bandpass)

% Identifying scans and on-axis clicks
max_ICI    = 500 ;          % Max ICI (ms) allowed to belong in a scan
min_clicks = 6 ;            % Minimum number of clicks required to include scan in analysis

% Graphics for click visualization tool (needs to reflect number of hydrophones)
symb = {'x';'o';'o';'x'};
siz  = [4 6 6 4];                 

if ~(length(channels)==length(channel_depths))
    disp('Warning: Reconfigure settings. Vector of channels and channel depths must be same size')
elseif ~(length(channels)==length(clip_level))
    disp('Warning: Reconfigure settings. Must contain a final sensitivity for each channel')
elseif ~(length(channels)==length(symb))
    disp('Warning: Reconfigure settings. Must contain a symbol for each channel')
elseif ~(length(channels)==length(siz))
    disp('Warning: Reconfigure settings. Must contain a marker size for each channel')
end

S.setupname  = setupname ;
S.datadir    = datadir ;
S.clickprefix= clickprefix ;
S.clip_level = clip_level ;
S.gain       = gain ;
S.soundspeed = soundspeed ;
S.temperature= temperature ;
S.channels   = channels ;
S.CH         = channel_center ;
S.R          = channel_depths(:) ;
S.max_ICI    = max_ICI ;
S.min_clicks = min_clicks ;
S.detect_filt= detect_filt ;
S.detect_bl  = detect_bl ;
S.clickdur   = detect_clickdur ;
S.localize_method = localize_method ;
S.NFFT       = NFFT ;
S.FFT_interp = FFT_interp ;
S.centroid   = centroid ;
S.filt       = filt ;
S.symb       = symb ;
S.siz        = siz ;