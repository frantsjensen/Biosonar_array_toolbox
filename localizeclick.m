function [LOC , RLpp , Tc , SIGNAL] = localizeclick(file,click,SETTINGS,userinput) ;

% [LOC,RLpp,T,SIGNAL] = localizeclick(file,click_sample,SETTINGS,userinput)
%
% Finds the source location of a click recorded on a vertical array
% using time-of-arrival differences (TOAD)
%
% Input:
%   file (path and file name)
%   click_sample (approximate sample of click)
%   SETTINGS is a structure with fields:
%       channels:           Channels used for source localization
%       R:                  Vector of depth of each receiving channel
%       CH:                 Central receiver (used by click detector)
%       soundspeed:         Sound speed in area
%       localize_method:    Time-of-arrival difference localization method
%                           Options: 'peak' (peak TOAD),'6db' (-6 dB shoulder),
%                           'correlation' (cross-correlation TOAD)
%       clickdur:           Expected click duration (in seconds)
%       detect_filt:        High-pass or bandpass filter (in Hz)
%       NFFT:               FFT size of analysis (for extracting click waveforms)
%       clip_level:         Clip level of all hydrophones used in click localization
%
% Output:
%   LOC: Source coordinate vector (m) in respect to vertical array and surface 
%   (x = horizontal distance, y = depth)
%   RLpp: Peak-to-peak sound pressure level (dB re 1 muPa)
%   T: Time of the click on all channels (samples since start of file)
%   SIGNAL: Matrix of signal waveform from all channels
%
% Requires multi-channel wav-file (compatible with wavexread) as input
%
% F. H. Jensen, 2013 (frants.jensen@gmail.com)

if nargin < 4,
    userinput = 0; % Default is no user input, fully automated
end

if nargin < 3,
    error(' Click localization needs settings input ')
end

channels    = SETTINGS.channels ;   % Channels used for click localization
R           = SETTINGS.R ;          % Depth of hydrophones used for click localization
CH          = SETTINGS.CH;          % Central hydrophone
c0          = SETTINGS.soundspeed ; % Sound speed
method      = SETTINGS.localize_method ; % Method {peak,correlation,6db}
siglength   = 3*SETTINGS.clickdur ;
HP          = SETTINGS.detect_filt ;
NFFT        = SETTINGS.NFFT ;
clip_level  = SETTINGS.clip_level ;

% Find file size and sample rate
[temp , fs] = wavexread ( file , [1 2] );   % Find fs

% Estimate length of window to search (based on size of array)
% maximum distance from click detection channel to rest of channels
% is calculated based on depth of each channel
fullwin       = floor(0.95*max(abs(R-R(CH)))*fs/c0); 

% Design filter
if length(HP)==2
    [BBB,AAA] = butter(2,HP/(fs/2)) ;               % Bandpass filter
else
    [BBB,AAA] = butter(2,HP/(fs/2),'high');         % Highpass filter
end

% Read and process the full click from four channels
sig = wavexread( file,click+fullwin*[-1 1]);

% Keep only channels used for localization
sig = sig(:,channels);

% Go through each channel and filter data
for i=1:length(R),
    sig(:,i) = sig(:,i) - mean (sig(:,i));      % Correct for DC
    sig(:,i) = filtfilt ( BBB,AAA,sig(:,i) );     % Filter Signal
    sig(:,i) = sig(:,i) .* tukeywin(length(sig),NFFT/length(sig));      % Apply Envelope
end

T_peak = zeros(1,length(R)) ;
T_6db  = zeros(1,length(R)) ;

% Find peak and -6 db shoulder of clicks
for i=1:length(R),
    hill(:,i)     = abs(hilbert(sig(:,i)));             % Peak detection
    [peaks,locs]  = findpeaks(hill(:,i),'minpeakheight',0.5*max(hill(:,i)),'sortstr','descend','minpeakdistance',10);
    if isempty(locs)
        T_peak(i) = find(hill(:,i)==max(hill(:,i)));
    elseif length(locs)==1,
        T_peak(i) = locs;
    else
        T_peak(i) = min(locs(1:2));
    end
    T_6db(i)      = min(find(hill(:,i)>0.5*max(hill(:,i))));
end

% Cross-correlate
model = sig(T_peak(CH)+[-round(siglength*fs):round(siglength*fs)],CH);
model = model.*tukeywin(length(model),0.3);
for i=1:length(R),
    x_corr          = xcorr ( sig(:,i),model);
    x_hill (:,i)    = abs ( hilbert(x_corr));
end

for i=1:length(R),
    [peaks,locs]    = findpeaks(x_hill(:,i),'minpeakheight',0.5*max(x_hill(:,i)),'sortstr','descend','minpeakdistance',10);
    if isempty(locs)
        xcpeak(i)  = find(x_hill(:,i)==max(x_hill(:,i)));
    elseif length(locs)==1,
        xcpeak(i)  = locs;
    else
        xcpeak(i)  = min(locs(1:2));
    end
end

T_corr = T_peak(CH) + xcpeak - xcpeak(CH);  % Cross-correlated TOADs in samples

% Select the proper time or arrival
switch method
    case 'peak' 
        T = T_peak ;
    case '6db'
        T = T_6db ;
    case 'correlation'
        T = T_corr ;
    otherwise
        T = T_corr ;
end

cont = 0 ;

while ~cont
    
    % Localize click using linear least squares
    t = ( T(2:end)' - T(1) )./fs;            % Convert TOADs to seconds relative to top hydrophone
    r = (R(2:end)-R(1));                     % normalized receiver coordinates
    A = 2 * [r t*c0^2];                      % source location A matrix, defined in Wahlberg et al. 2001 and this ms
    B = - (t.^2)*(c0^2) + r.^2;              % source location b column vector, defined in Wahlberg et al. 2001 and this ms
    m = A\B;                                 % source solution solved by least squares as in Wahlberg et al. 2001 and this ms
    so(2) = m(1)+R(1);                       % source depth coordinate
    so(1) = sqrt((c0*m(2))^2 - m(1)^2 );     % source horizontal coordinate

    % Rename location
    LOC = real(so) ;

    % Calculate range to central hydrophone
    % range = ( sqrt( so(1)^2 + (so(2)-R(CH)).^2) )';

    % Correct time of arrivals back to file samples
    Tc = round(T) + click-fullwin-1 ;

    % Find short window of click on each channel
    SIGNAL = zeros(NFFT,length(R));
    for i=1:length(R),
        shortsig = wavexread( file,Tc(i)+[-NFFT/2+1 NFFT/2]) ;
        SIGNAL(:,i) = shortsig(:,i)-mean(shortsig(:,i)) ;
    end

    % Find peak-to-peak RL for each channel
    for i=1:length(R),
        RLpp(i) = 20*log10(max(SIGNAL(:,i))-min(SIGNAL(:,i)))+clip_level(i);
    end
    
    % Check if user input is required. If it is, plot localization
    if userinput

        % Define colors for plot - only 8 colors defined per default
        if length(R)>8
            error('Update channelcolors.m with more color schemes to allow for more than 8 channels')
        else
            CMAT = channelcolors(length(R));
        end

        %%%%%%%%%%%%%%%           START PLOT            %%%%%%%%%%%%%%%%
        figure(11), set(gcf,'Name','Sound source localization'), clf, subplot(3,2,1),hold on, box on
        for i = 1:length(R), plot( [1:length(sig)] , sig(:,i) ,'k','color',CMAT(i,:),'LineWidth',2), legg{i} = ['Ch. ' num2str(i)]; end
        axis ( [min(T_peak)-round(siglength*fs) max(T_peak)+3*round(siglength*fs) 0.95*min(min(sig)) 1.05*max((max(sig)))])
        hleg = legend(legg,'Location','NorthEast'); 
        title('Signal Waveforms');
        set(hleg,'box','off')
        xlabel('Time, samples')
        ylabel('Amplitude')

        subplot(3,2,2), hold on, box on,
        for i = 1:length(R), 
            h(i)=plot( [1:length(hill(:,i))] , hill(:,i) ,'k','color',CMAT(i,:),'LineWidth',2);
            plot( [T_peak(i) T_peak(i)],[0 max(max(hill(:,i)))],'k--','color',CMAT(i,:),'LineWidth',2)
        end
        axis ( [min(T_peak)-round(siglength*fs) max(T_peak)+3*round(siglength*fs) 0 max(max(hill))])
        hleg = legend(h,legg,'Location','NorthEast'); title('Signal Envelopes');
        set(hleg,'box','off')
        xlabel('Time, samples')
        ylabel('Amplitude')

        subplot(3,2,3), hold on, box on
        plot ( [1:length(sig(:,CH))],sig(:,CH)/max(sig(:,CH)),'k',T_peak(CH)-round(siglength*fs)+[0:length(model)-1],tukeywin(length(model),0.3),'r','LineWidth',2), 
        axis ( [T_peak(CH)+[-1 1]*round(siglength*fs) 1.05*min(sig(:,CH))/max(sig(:,CH)) 1.05] )
        title( ['Cross Correlation model (Ch ' num2str(CH) ') and Tukey window (red)'])
        xlabel('Time, samples')
        ylabel('Normalized amplitude')

        if ~strcmp(method,'peak') & ~strcmp(method,'6db'),
            subplot(3,2,4), hold on, box on
            for i = 1:length(R), 
                h(i)=plot( [1:length(x_hill(:,i))],x_hill(:,i),'k','color',CMAT(i,:),'LineWidth',2);
                plot(T_corr(i)+(xcpeak(CH)-T_peak(CH))*ones(2,1),[0 max(x_hill(:,i))],'k--','color',CMAT(i,:),'LineWidth',2), 
            end    
            axis ( [min(xcpeak)-round(siglength*fs) max(xcpeak)+3*round(siglength*fs) 0 1.05*max(max(x_hill))])
            xlabel('Cross correlation samples')
            ylabel('Correlation amplitude')
        end
        hleg = legend(h,legg,'Location','NorthEast');  title('Cross Correlation Envelopes');
        set(hleg,'box','off')

        % Localization plot
        subplot(3,2,5), hold on, box on
        step = 0.01 ;
        ax = [-1 round(so(1)+5) min([R(4) so(2)]-2) 1] ;
        for i=1:length(r),                       
            a = c0*t(i)/2;                                          % the three parameters a, b, and c used to define the hyperbola curve in line 29
            c = r(i) / 2;                
            b = sqrt(c^2 - a^2);
            if imag(b) == 0,                                        % check if TOAD render physical hyperbola
                if t(i) ~= 0                     
                    y = -sign(r(i))*sign(t(i))*[step:step:2*(ax(4)-ax(3))] ;  % equally spaced vector for which to calcualate a hyperbola below
                    x = b * sqrt( y.^2 / a^2 - 1);                  % equation of a hyperbola curve with the parameters a and b defined above
                else
                    y = zeros(1,2);                                 % plot straight line instead if TOAD is zero
                    x = [0 ax(2)];
                end
                ind = min(find(imag(x)==0));                        % find index in x vector where hyperbola starts
                x = [ fliplr(-x(ind:end)) x(ind:end) ]';            % extend hyperbola to both positive and negative x values
                y = [ fliplr(y(ind:end))  y(ind:end)]' + c + R(1);  % extend y vector to two-sided hyperbola curve
                plot(x,y,'k','color',CMAT(i,:))                     % plot hyperbola
            end
        end
        plot([ax(1):0.01:ax(2)],0.1*sin(2*pi*[ax(1):0.01:ax(2)]),'b','LineWidth',2)
        plot(so(1),so(2),'r*'); plot([0 0],[min(R) 0],'k','LineWidth',2); 
        plot(zeros(length(R),1),R,'ko','MarkerSize',4,'MarkerFaceColor','k')
        axis(ax)
        title('TOAD Localization')
        xlabel('Horizontal distance, m')
        ylabel('Depth, m')
        
        subplot(3,2,6), hold on, box on
        for i = 1:length(R), plot(1000*[1:length(SIGNAL)]/fs,SIGNAL(:,i),'k','color',CMAT(i,:),'LineWidth',2), legg{i} = ['Ch. ' num2str(i)]; end
        axis ( [0 , 1000*NFFT/fs , 1.1*min(min(SIGNAL)) , 1.1*max((max(SIGNAL)))])
        hleg = legend(legg,'Location','SouthEast'); title('Synchronized clicks');
        set(hleg,'box','off')
        xlabel('Time, ms'), ylabel('Rel. Amplitude')

        % Change font name and font size
        fontname = 'helvetica' ;
        fontsize = 12 ;
        
        % Find handles for all axes and labels within current figure
        handle_axes = get(gcf,'children') ;
        handle_labels = get(handle_axes,{'XLabel','YLabel'});
        handle_titles = get(handle_axes,{'Title'});

        % Change font on those axes and labels
        set(handle_axes,'FontName',fontname,'FontSize',fontsize-2)
        set([handle_labels{:,:}],'FontName',fontname,'FontSize',fontsize)
        for j=1:length(handle_titles), 
            set(handle_titles{j},'FontName',fontname,'FontSize',fontsize+2), 
        end
        
        % Enter user input choices
        userchoice = centermenu(' Accept localization or change method ','Accept','Cross correlation','Peak envelope','-6 dB start of envelope');
        switch userchoice
            case 2
                disp('Cross correlation TOADs')
                T = T_corr ;
            case 3
                disp('Peak envelope TOADs')
                T = T_peak ;
            case 4
                disp('-6 dB start of click TOADs')
                T = T_6db ;
            otherwise
                cont = 1;
        end     
        
    else % If user input not required, end while loop
        cont = 1 ;
    end
    
end