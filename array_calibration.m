% Click localization and precision estimator
% Newest revision 10.08.2013 by Frants Jensen
% Used for calculating source localization precision and accuracy
%
% Adjust variables manually in lines 20 to 39.
%
% Progress of Script:
% A: Calibration
% B: Progressive sound file analysis - includes:
% b1: Click Detection
% b2: Click Ranging
% -Peak of envelope or cross correlation choice (cross correlation best)
% C: Analysis, quantification of performance
%
% F. H. Jensen 2013 (frants.jensen@gmail.com)

%%%%%%%%%%%%%%%%    VARIABLES    %%%%%%%%%%%%%%%%%%%
clear, close all,

library         = 'C:\Users\FHJ\Dropbox\SousaClick\';
filename        = 'Reson_10m_40dB_3_1.0V(peak)__04_02_2014_17_04_02.wav';
file            = [library filename];

threshold       = 0.012;    % Define click threshold - change this depending on level of received calibrations
channel         = 5;        % Define channel to be searched for clicks (center of array)
jump            = 14000;    % Define number of samples to pass after click detected
                            % -1500 samples equal 1875 usec (or click freq 400 click/sec)
                            % You can define this based on known repetition
                            % rate of calibration signal
N1              = 1;        % Define first sample to be analysed
HP_detect       = 50000;    % HP filter for detection, in Hz 
HP_localize     = 10000;    % HP filter for localization, in Hz
chunk           = 20000;    % Number of samples to be read at once
clip_level      = 190;      % Define peak clip level [220+14 dB re V - AMP dB ] = 164 with 70dB Gain, 174 with 60dB gain and 184 with 50dB gain
c0              = 1434;     % Speed of sound [m/s]
Temp            = 1;        % Av. temperature from outer harbour
cross_correlate = 1;        % If true (1) uses cross-correlated TOADs instead of Peak-detected
R = [-1 -1.75 -2.5 -3.25 -4]' ; % List of receiver depth coordinates, must be column, all negative
receivers = [1 3 5 6 7];        % receivers to use for localization

%%%%%%%%%%%%%%%% b1: Click Detection    %%%%%%%%%%%%%%%%%%%
n_clicks = 0;

[file_size]         = wavexread ( file , 'size' );% Find total size of sound
[sig_chunk , fs]    = wavexread ( file , [1 2] ); % Find fs
if file_size(2)~=6, 
    disp(['Warning: ' num2str(file_size(2)) ' channels!']), 
end
file_size = file_size(1);                           % Discard number of channels

[BBB,AAA] = butter ( 4 , HP_detect / (fs/2) , 'high');  % Define filter
segment_length = chunk;

tic
%The following while-loop will run until the end of the file
while N1 < file_size - jump

    %Define part of sound file to read
    if N1+chunk-1 > file_size,
        [sig_chunk,fs] = wavexread (file , [N1 file_size]);
        segment_length = file_size-N1;
    else
        [sig_chunk,fs] = wavexread (file , [N1 N1+chunk-1]);
    end
    
    % Clear memory of unnecessary data
    sig_chunk = sig_chunk(:,channel); 

    % Warn if data are clipped
    if max(abs(sig_chunk))>=0.99
        disp ('warning - input data might be clipped' );
        disp( ['progress: N1 = ' num2str(N1) ] )
    end

    int = 1;    
    sig_chunk = sig_chunk - mean (sig_chunk);
    sig_chunk = filter (BBB , AAA , sig_chunk);

    % Look for clicks through chunk
    while int < segment_length,

        % Click located if absolute value exceeds threshold
        if  abs ( sig_chunk (int) ) > threshold,

            % if click is not too close to last click, save position of click
            if ~exist('click')
                n_clicks = n_clicks + 1;
                click(n_clicks) = N1 + int; % Varies w. about 20 samples from exact time                
            elseif (N1 + int > click(end) + jump)
                n_clicks = n_clicks + 1;
                click(n_clicks) = N1 + int; % Varies w. about 20 samples from exact time
            end
 
            % Once peak has been detected, skip a few samples 
            % before continuing click detector:
            int = int + jump;
        end
        int = int + 1; % Can be set to int+2???
    end

    % Once end of chunk has been reached, proceed with next chunk:
    N1  = N1 + chunk;
end
toc

if ~exist('click'), error('No clicks detected in file - adjust threshold?'),end
disp(['Click Detection Complete for file ' file])
disp([num2str(max([N1-chunk click(end)])) ' samples out of ' num2str(file_size) ' analysed and ' num2str(n_clicks) ' clicks detected.']) 

save( 'calibration_test.mat' , 'click' , 'n_clicks' , 'fs' , 'file' , 'Temp' , 'c0' , 'clip_level' , 'cross_correlate','HP_localize','R','receivers');
disp(' Saved calibration results to calibration_test.mat')
disp(' Rename and save to data directory ')

%%%%%%%%%%%%%%%% b2: Click Analysis     %%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%Parameters%%%%%%%%%
clear
load('calibration_test.mat')

full_win        = ceil(fs*3/1500);  % Analyses 2 ms (3m/1500[m/s]) on either side of click detected
fft_win         = 128;

% Defining filter
[BBB,AAA] = butter( 2 , HP_localize/(fs/2) , 'high');

% Configure ranging plot
pattern= ['b-';'g-';'r-';'c-'; 'm-' ; 'y-'];

if click(1)<full_win
    click(1)=[];
    n_clicks=n_clicks-1;
end

if n_clicks > 1000, n_clicks = 1000; end

for click_j = 1:n_clicks;
    sigtemp = wavexread( file, click(click_j) + full_win * [-1 1] );
   
    sig=zeros(length(sigtemp),length(receivers));
    for i = 1:length(receivers)
        sig(:,i) = sigtemp(:,receivers(i));
    end
    
    clear hill T_peak x_corr hill T_xcorr kryds_max
    for i=1:length(R),
        sig(:,i) = sig(:,i) - mean (sig(:,i)) ;
        sig(:,i) = filter ( BBB,AAA,sig(:,i) );     % Check filter magnitude response
        sig(:,i) = sig(:,i).*tukeywin(length(sig(:,i)),128/length(sig(:,i)));
    end
    
    %Ignore potential surface reflection - only applies to artificial porpoise clicks
    clear hill T_* temp*
    for i=1:length(R),
        hill (:,i)    = abs  ( hilbert(sig(:,i)) );      % Peak detection
        [peaks,locs]  = findpeaks(hill(:,i),'minpeakheight',3*std(hill(:,i)),'sortstr','descend','minpeakdistance',10);
        if isempty(locs)
            T_peak(i) = find(hill(:,i)==max(hill(:,i)));
        elseif length(locs)==1,
            T_peak(i) = locs;
        else
            T_peak(i) = min(locs(1:2));
        end
    end

	% Cross-correlate with direct click only - window off surface reflections
	han_w = hanning (64); han_w = han_w(1:32);
	han_w2 = [han_w ; ones(16,1) ; flipud(han_w)];
	
	if T_peak(3)+30 > length(sig),
        temp_sig = sig(:,3);
        % disp ('warning, may have to choose bigger window')
	else
        cutoff = zeros(length(sig),1); 
        cutoff([-39:40]+T_peak(3)) = han_w2;
        temp_sig = sig(:,3).*cutoff;
	end
    
    clear x_* 
    for i=1:length(R),
        x_corr        = xcorr ( sig(:,i),temp_sig);
        x_hill (:,i)  = abs ( hilbert(x_corr));
        [peaks,locs]  = findpeaks(x_hill(:,i),'minpeakheight',3*std(x_hill(:,i)),'sortstr','descend','minpeakdistance',10);
        if isempty(locs)
            T_xcorr(i) = find(x_hill(:,i)==max(x_hill(:,i)));
        elseif length(locs)==1,
            T_xcorr(i) = locs;
        else
            T_xcorr(i) = min(locs(1:2));
        end
    end
    
    if click_j > 1,
        IPI(click_j)  = 1000*( click(click_j) - click(click_j-1) )/fs;
    end
       
    % Configure hydrophone depths and time differences
    T = T_peak(1) + T_xcorr - T_xcorr(1);    % Cross-correlated TOADs in seconds
    %T = T(receivers);
    %R = R(receivers);
    
    % Localize click    
    t = ( T(2:end)' - T(1) )./fs;            % Peak-detected TOADs in seconds
	r = (R(2:end)-R(1));                     % normalized receiver coordinates
	step = .01;                              % step size in hyperbola plots
	A = 2 * [r t*c0^2];                      % source location A matrix, defined in Wahlberg et al. 2001 and this ms
	B = - (t.^2)*(c0^2) + r.^2;              % source location b column vector, defined in Wahlberg et al. 2001 and this ms
    m = A\B;                                 % source solution solved by least squares as in Wahlberg et al. 2001 and this ms
	so(2) = m(1)+R(1);                       % source depth coordinate
	so(1) = sqrt((c0*m(2))^2 - m(1)^2 );     % source horizontal coordinate
    
    % Used for box plot of localizations - axes control
    % Retain the biggest size of the localization plot    
    if ~exist('ax')                          
        ax = [-1 round(so(1)+5) min([R(4) so(2)]-2) 1];
    else
        ax = [-1 max([round(so(1)+5) ax(2)]) min([R(4)-2 so(2)-2 ax(3)]) 1];
    end
    
    % Range between source and each hydrophone is calculated
    range = ( sqrt( so(1)^2 + (so(2)-R).^2) )';  
    
    % Add these coords to the localization plot

    %%%%%%%%%%%%%%%START PLOT%%%%%%%%%%%%%%%%
    if any (click_j == plot_click)
        subplot(2,2,3), hold on, 
        for i = 1:size(x_hill,2), plot( [1:length(x_hill(:,i))] , x_hill(:,i) , pattern(i,:) , [T_xcorr(i) T_xcorr(i)],[0 max(max(x_hill))],pattern(i,:) ),end
        axis ( [min(T_xcorr)-ceil(0.001*fs) max(T_xcorr)+ceil(0.001*fs) 0 max(max(x_hill))])
        title('Cross-correlation Envelopes and detected peaks - b g r c m y')

        subplot(2,2,4),hold on
        plot(so(1),so(2),'r*');                  % plot source coordinates
        % Add hyperbolas to the localization plot (currently deactivated)
		for i=1:length(r),                       
            a = c0*t(i)/2;                       % the three parameters a, b, and c used to define the hyperbola curve in line 29
            c = r(i) / 2;                
            b = sqrt(c^2 - a^2);
            if imag(b) == 0,                     % check if TOAD render physical hyperbola
                if t(i) ~= 0                     
                    y = -sign(r(i))*sign(t(i))*[step:step:2*(ax(4)-ax(3))] ;  % equally spaced vector for which to calcualate a hyperbola below
                    x = b * sqrt( y.^2 / a^2 - 1); % equation of a hyperbola curve with the parameters a and b defined above
                else
                    y = zeros(1,2);             % plot straight line instead if TOAD is zero
                    x = [0 ax(2)];
                end
                ind = min(find(imag(x)==0));    % find index in x vector where hyperbola starts
                x = [ fliplr(-x(ind:end)) x(ind:end) ]'; % extend hyperbola to both positive and negative x values
                y = [ fliplr(y(ind:end))  y(ind:end)]' + c + R(1); % extend y vector to two-sided hyperbola curve
                plot(x,y)                       % plot hyperbola
            end
		end
        for i=1:length(r)-1, plot ( loc_xy (i,1), loc_xy (i,2), 'bo'), end
        axis ([-1 round(so(1)+5) min([R(end) so(2)]-2) 1]);
        title('Localizations')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
    SOURCE_XY (click_j,:) = so;
    RANGE     (click_j,:) = range;
end

SOURCE_X = SOURCE_XY(:,1) ;
SOURCE_Y = SOURCE_XY(:,2) ;

h2=figure(2); clf,
set(h2,'Color',[1 1 1]);
set(h2,'PaperUnits','normalized')
set(h2,'Position',[10 50 1000 600])
subplot(2,1,1),box on, 
plot( real(SOURCE_XY(:,1)) , real(SOURCE_XY(:,2)) , 'rx' ), 
set(gca, 'FontName','helvetica','FontSize',12)
title ('Source coordinates','FontName','helvetica','FontSize',14), 
xlabel('Horizontal range, m','FontName','helvetica','FontSize',14), 
ylabel('Estimated depth, m','FontName','helvetica','FontSize',14)

subplot(2,1,2), box on
plot( [1:length(SOURCE_XY)] , real(SOURCE_XY(:,1) )), 
set(gca, 'FontName','helvetica','FontSize',12)
title ('Drift','FontName','helvetica','FontSize',14), 
xlabel('Click number (i.e. time)','FontName','helvetica','FontSize',14), 
ylabel('Range from array','FontName','helvetica','FontSize',14)

if any(any(imag(SOURCE_XY)))
	temp = find(~imag(SOURCE_X));
    disp ([' Filtered ' num2str(length(SOURCE_X)-length(temp)) ' from data due to imaginary parts'])
    SOURCE_X = SOURCE_X(temp) ;
    SOURCE_Y = SOURCE_Y(temp) ; 
end

min_range = input ('Filter for range? Enter minimum range (0 = no filter): ');
if min_range 
	temp = find(abs(SOURCE_X)>min_range);
	disp([' Filtered ' num2str(length(SOURCE_X)-length(temp)) ' detected clicks from data (min range)'])
    SOURCE_X = SOURCE_X(temp);
	SOURCE_Y = SOURCE_Y(temp);
end

max_range = input ('Filter for range? Enter maximum range (0 = no filter): ');
if min_range 
	temp = find(abs(SOURCE_X)<max_range);
	disp([' Filtered ' num2str(length(SOURCE_X)-length(temp)) ' detected clicks from data (max range)'])
	SOURCE_X = SOURCE_X(temp);
	SOURCE_Y = SOURCE_Y(temp);
end

range = input ('What is the real range? ');
depth = -2;


%%%%%%%%%%%%%%% PRESENT DATA %%%%%%%%%%%%%%%
h2=figure(2); clf,
set(h2,'Color',[1 1 1]);
set(h2,'PaperUnits','normalized')
set(h2,'Position',[10 50 1000 600])
subplot(2,1,1), box on
plot( SOURCE_X , SOURCE_Y , 'rx' ), 
set(gca, 'FontName','helvetica','FontSize',12)
title ('Source coordinates','FontName','helvetica','FontSize',14), 
xlabel('Horizontal range, m','FontName','helvetica','FontSize',14), 
ylabel('Estimated depth, m','FontName','helvetica','FontSize',14)

subplot(2,1,2), box on
plot( [1:length(SOURCE_X)] , SOURCE_X ,'kx'), 
set(gca, 'FontName','helvetica','FontSize',12)
title ('Drift','FontName','helvetica','FontSize',14), 
xlabel('Click number (i.e. time)','FontName','helvetica','FontSize',14), 
ylabel('Range from array','FontName','helvetica','FontSize',14)

disp(['N = ' num2str(length(SOURCE_X))]);

RMSE = sqrt( sum((SOURCE_X-range).^2)/length(SOURCE_X) ) ;

%%%%%%%%%%%%%%% Least-Square-derived source positions %%%%%%%%%%%%
disp(' ')
disp([' CALIBRATION RESULTS AT ' num2str(range) 'm RANGE:'])
disp(' ')
disp([' Number of detected clicks (N): ' num2str(length(SOURCE_X))])
disp([' Horisontal Mean range = ' num2str(mean(SOURCE_X)) ])
disp([' Horizontal Standard deviation = ' num2str(std(SOURCE_X))])
disp([' Horisontal RMS error, m = ' num2str(RMSE)])
disp([' Horisontal RMS error, % = ' num2str(100*RMSE/range)])
