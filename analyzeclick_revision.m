function [] = analyzeclick_revision(file,Tc,LOC,ICI,SETTINGS,currentfile)
% [] = analyzeclick_revision(file,Tc,LOC,ICI,SETTINGS,currentfile)
%
% Parent script reviseanalysis.m: Finds the source parameters of biosonar 
% clicks that have already been extracted from sound files, using
% information about source files and click locations contained in
% previously saved *params.mat files.
%
% The parent script, reviseanalysis.m, will call analyzeclick_revision.m
% for each existing *params.mat file. Thus, these two scripts will allow you
% to change analysis parameters (such as FFT size) or add analysis
% components (such as different metrics) after clicks have already been
% extracted, and will update click parameter files to reflect the new
% analysis settings.
%
% All parameters extracted following equations in Madsen and Wahlberg 2007
% and W. W. L. Au (1993): Sonar of Dolphins
%
% F. H. Jensen, 2013 (frants.jensen@gmail.com)


% Find variables from settings file
R           = SETTINGS.R ;
NFFT        = SETTINGS.NFFT ;           % Desired FFT size
fft_int     = SETTINGS.FFT_interp ;     % FFT interpolation
HP          = SETTINGS.filt ;           % Filter settings, Hz
centroid    = SETTINGS.centroid ;       % Centroid frequency in kHz - for calculating absorbtion
tempr       = SETTINGS.temperature ;    % Temperature of water, for calculating absorbtion.
clip_level  = SETTINGS.clip_level ;     % Clip level of all transducers
datadir     = SETTINGS.datadir ;        % Where to store output
clickprefix = SETTINGS.clickprefix ;    % Prefix for output files

hydcorr = dlmread('hydrophone_correction.txt');

% Find short window of click on each channel
sig = zeros(NFFT,length(R));
for j=1:length(R),
    % To correct hydrophone frequency response, we load a slightly larger
    % window, implement correction, and re-select a window corresponding to
    % the nominal FFT size
    [sigtemp fs] = wavexread(file,Tc(j)+[-64+1 64]) ;
    sigtemp = ifft(fft(sigtemp(:,j))./hydcorr(:,j)) ;
    sig(:,j) = sigtemp(64+[-NFFT/2+1:NFFT/2]) ;
    
    % Default code:
    %    [sigtemp fs] = wavexread(file,Tc(j)+[-NFFT/2+1 NFFT/2]) ;
    %    sig(:,j) = sigtemp(:,j) ;
    
end

% Check if waveforms are clipped
if max(max(abs(sig)))>=0.95
    disp ('warning - input data likely clipped' );
end

% Calculate absorption loss based on temperature and expected centroid frequency
A_abs = 48.83e-8 + 65.34e-10*tempr;
B_abs = 1.55e7 * (tempr+273.1) * exp(-3052/(tempr + 273.1));
alpha = A_abs*B_abs*(1000*centroid)^2 / (B_abs^2 + (1000*centroid)^2); %Unit dB/m

% Design filter
if length(HP)==2
    [BBB,AAA] = butter(4,HP/(fs/2)) ;               % Bandpass filter
else
    [BBB,AAA] = butter(4,HP/(fs/2),'high');         % Highpass filter
end

% Go through channels and filter signals
for j=1:length(R),
    sig(:,j) = sig(:,j)-mean(sig(:,j));             % DC correction
    sig(:,j) = filter ( BBB,AAA,sig(:,j) );       % Filter
    sig(:,j) = sig(:,j) .* hanning(length(sig));       % HANNING WINDOW
    % Default code:
    % sig(:,j) = sig(:,j) .* hanning(length(sig));       % HANNING WINDOW
    % Possibly correction for 200 Hz LP filter in amplifier box
end

% Zeropad for speed:
fc = zeros(1,length(R)); fp = zeros(1,length(R)); bw3db = zeros(1,length(R)); 
bw10db = zeros(1,length(R)); rmsbw = zeros(1,length(R)); xpsd = zeros(0.5*NFFT*fft_int,length(R));

% Estimate spectral parameters (following Madsen and Wahlberg 2007)
for j=1:length(R),
    fsig     = abs(fft(sig(:,j),NFFT*fft_int));                 % High-quality Sinc fft interpolation
    fsig     = fsig(1:length(fsig)/2);                          % Drop second half of amplitude spectrum
    xpsd(:,j)= 2 * (1/fs/(NFFT*fft_int).*abs(fsig).^2) ;        % Power spectral density estimate
    fsig     = fsig/max(fsig);                                  % Normalize amplitude spectrum
    f        = [0:length(fsig)-1]'*fs/(2*10^3*(length(fsig)));  % Calculate frequency vector
    fc(j)= sum(f.*fsig.^2)/sum(fsig.^2);                        % Centroid frequency, kHz
    rmsbw (j)= sqrt(sum((f-fc(j)).^2.*fsig.^2 ) / sum(fsig.^2));% RMS centralized bandwidth, kHz
    fp(j)    = f(find(fsig==max(fsig)));                        % Peak frequency in kHz
    f_index  = find(fsig.^2>10^(-3/10));                        % Find index of power spectrum within -3dB of peak power
    bw3db(j) = (f(max(f_index))-f(min(f_index)));               % -3 dB BW in kHz
    f_index  = find(fsig.^2>10^(-10/10));                       % Find index of power spectrum within -10dB of peak power
    bw10db(j)= f(max(f_index))-f(min(f_index));                 % -10 dB BW in kHz
    
    %Constructing time axis
    time_axis = 10^6 * ([1:length(sig)]) / fs;
end

% Define colors for plot - only 8 colors defined per default
if length(R)>8
    error('Update channelcolors.m with more color schemes to allow for more than 8 channels')
else
    CMAT = channelcolors(length(R));
end

figure(12), clf, set(gcf,'Name','Click parameter analysis')
peakpsd = max(10*log10(max(xpsd))+clip_level) ;
for j=1:length(R)
     % Plot signals and frequency spectrum
    subplot(2,1,1), hold on, box on
    plot (time_axis,sig(:,j),'k','LineWidth',2,'color',CMAT(j,:)); legg{j} = ['Ch. ' num2str(j)];
    axis ( [0 1e6*NFFT/fs min(min(sig)) max(max(sig))] ),xlabel('usec'), ylabel('Rel. amplitude')
    
    subplot(2,1,2), hold on, box on
    h(j)=plot (f,10*log10(xpsd(:,j))+clip_level(j),'k','Linewidth',2,'color',CMAT(j,:));
    plot ([-0.5 0.5]*rmsbw(j)+fc(j), [1 1]*max(10.*log10(xpsd(:,j))+clip_level(j))-20, 'r' ,'Linewidth',3,'Color',CMAT(j,:)),
    axis ( [0 max(f) [-40 0]+10*ceil(0.1*(peakpsd))] ), xlabel('kHz'), ylabel('PSD, dB re 1muPa2/Hz')
end
subplot(2,1,2), hleg = legend(h,legg,'location','NorthEast'); set(hleg,'box','off')

% Find exact range and angle of arrival to hydrophones
range = localization2range(LOC,R);

% Calculate and store source parameters.
% ASLpp = RLpp + TL, TL = 20*log(range) + alpha*range
for j=1:length(R)
    % -10 dB amplitude criterion
    [peak , rms , window, flux] = RL ( sig(:,j), fs, clip_level(j) , 1);
    RLpp       (j) = peak;
    ASLpp      (j) = peak + 20*log10( range(j) ) + alpha*range(j);
    RLrms_amp  (j) = rms;
    ASLrms_amp (j) = rms  + 20*log10( range(j) ) + alpha*range(j);
    ASLflux_amp(j) = flux + 20*log10( range(j) ) + alpha*range(j);
    WINDOW_amp (j) = window;

    % 95% Energy Criterion
    [peak , rms , window, flux] = RL ( sig(:,j), fs, clip_level(j) , 2);
    RLrms_en   (j) = rms;
    ASLrms_en  (j) = rms  + 20*log10( range(j) ) + alpha*range(j);
    ASLflux_en (j) = flux + 20*log10( range(j) ) + alpha*range(j);
    WINDOW_en  (j) = window;
end

% Find approximate angle of arrival to hydrophones 
% (assuming on-axis receiver is perfectly on-axis)
theta = localization2angle(LOC,R);

SIGNAL (:,:) = sig(:,:);
SPEC(:,:)    = xpsd(:,:);
FREQ         = f ;          % Frequency vector for power spectral density estimates

Fp           = fp(:);
Fc           = fc(:);
RMSBW        = rmsbw(:);
BW_3db       = bw3db(:);
BW_10db      = bw10db(:);

SOURCE_LOC   = LOC ;
RANGE        = range(:);
ANGLE        = theta(:);

ch=find(ASLpp==max(ASLpp));
fprintf(' Channel: \t\t\t%1.0f \n',(ch)) ;
fprintf(' ASLpp: \t\t\t%3.0f dB\n',ASLpp(ch)) ;
fprintf(' Fc: \t\t\t\t%2.0f kHz\n',Fc(ch)) ;
fprintf(' -3dB BW: \t\t\t%2.0f kHz\n',BW_3db(ch)) ;
fprintf(' Dur: \t\t\t\t%2.0f µs\n',WINDOW_amp(ch)) ;

save (currentfile, 'file','Tc','fs','ICI','SIGNAL','SPEC', 'FREQ','ch', ...
    'RLpp', 'ASLpp', 'RLrms_amp', 'ASLrms_amp', 'ASLflux_amp', 'WINDOW_amp', ...
    'RLrms_en', 'ASLrms_en' , 'ASLflux_en' , 'WINDOW_en' , ...
    'RMSBW', 'BW_3db' , 'BW_10db' , 'Fp' , 'Fc', 'RANGE', 'ANGLE', ...
    'SOURCE_LOC','SETTINGS')

disp(' ')
disp('Data sucessfully saved to: ') 
disp(currentfile)
disp(' ')