function [EPR]=estimateaperture(file,maxrange)
%
% EXPERIMENTAL!
%
% function []=estimateaperture(file,maxrange)
% Estimate the equivalent aperture and corresponding beam pattern for
% individual on-axis biosonar clicks using a parametric spectral fit
% (Jensen et al. submitted, JEB)
%
% Input:
% - file        file name (if file is in current directory) or 
%               path and file name of compiled click parameter file
%               (output from compileclickparams.m)
% - maxrange    Maximum range criterion for clicks (typically 20m)
%
% Output:
% - EPR         estimated equivalent piston radius for each click
% 
% Requires interpolate2.m

if nargin<2,
    maxrange=20;
end


% First, load compiled click parameter file
%D=load('C:\Matlab_Scripts\biosonar_array_toolbox\boto_Nov2014.mat') ;
D=load('file') ;

% Define variables
pistsize=2:0.01:20; % test these aperture diameters (cm)
ip = 8 ;            % Waveform interpolation factor for convolution


try 
    SETTINGS = toolbox_settings();
    c0 = SETTINGS.soundspeed;
catch
    c0 = 1524 ; % Default sound speed (used for convolution)
end

% Clip level difference (adjust these if click waveform not corrected)
clipleveldiff = 0+[0 0 0 0 0 0 0] ;

% Take the variables we need for this analysis
k=find(D.RANGE<maxrange);   % Index of signals produced within 20m confidence
ANGLES = D.ANGLES(k,:);     % Estimated angles
RELASL = D.RELASL(k,:);     % Back-calculated source levels
range = D.RANGE(k,:);       % Estimated range to on-axis hydrophone
Fc = D.Fc(k,:);             % Estimated on-axis centroid frequency
SIGNAL = D.SIGNAL(k,:);     % Waveform of signal recorded closest to acoustic axis
SIGS = D.SIGS(k);           % Waveform of all recorded signals
fs = D.fs(1) ;              % Sampling rate for waveforms
N = length(SIGNAL(1,:));    % FFT size for spectral analysis

% Start wait bar
h=waitbar(0,' Making parametric spectral fit for on-axis clicks');

% Go through each click
for j=1:size(SIGNAL,1),
    % For each click, find the waveform recorded closest to acoustic axis
    thissig = SIGNAL(j,:) ;
    
    % Since waveforms are not corrected for sensitivity difference,
    % we will correct them now
    k = find(RELASL(j,:)==max(RELASL(j,:))) ;
    thissig = thissig .*10.^(clipleveldiff(k)./20);
            
    for k=1:size(ANGLES,2),
        
        % For each channel, find the observed click waveform
        x = SIGS{j}(:,k) ;
        x = x(:).*rectwin(length(x)) .*10.^(clipleveldiff(k)./20);
                
        % Then calculate the observed amplitude spectrum for this receiver
        fx     = abs(fft(x,N));           % Calculate amplitude spectrum
        fx     = fx(1:length(fx)/2);      % Drop second half of amplitude spectrum
        freq     = (1/1000)*(fs/2)*[0:length(fx)-1]'/(length(fx)) ;
        
        %figure(1), hold on, plot(freq,10*log10(fx),'k','color',[k/size(ANGLES,2) 0 0])
        
        % Find the corrected angle of incidence for this receiver
        thisangle = ANGLES(j,k) ;
        
        % Parametric fit: 
        % Go through each piston size, model the expected amplitude spectrum 
        % at the recorded angle, and calculate error as the sum of squared 
        % deviations between expected and observed amplitude spectrum
        for pistindex = 1:length(pistsize)
            
            % Interpolate click
            fsint=ip*fs; s2 = interpolate2(thissig,ip) ;
            
            % Convolve with impulse response of circular piston
            xconv = conv(s2,imppist(abs(thisangle),pistsize(pistindex)/100,c0/1000,fsint/1000),'same');
            
            % Downsample to original sample rate
            xsim = xconv(ip:ip:end) ;
            
            % Calculate expected (simulated) amplitude spectrum
            fxsim     = abs(fft(xsim(:).*rectwin(length(xsim)),N));  % Calculate amplitude spectrum
            
            fxsim     = fxsim(1:length(fxsim)/2);        % Drop second half of amplitude spectrum
                     
            % Calculate total squared deviation (sum of squared error): sum(e-o)^2
            SSE(pistindex,k) = sum((fxsim-fx).^2) ; %SSE
        
        end       
    end
    
    % Find channels that are within 2-30 degrees off-axis
    kch = find(abs(ANGLES(j,:))>2 & abs(ANGLES(j,:))<30) ;

    % Calculate the total sum of squared errors as a function of piston size for those channels
    SSEtot = sum(SSE(:,kch),2) ;
       
    % Find the piston size that minimizes total squared error
    [test minepr] = min(SSEtot) ;
    EPR(j) = 0.5*pistsize(minepr) ;    
    
    % Calculate -3 dB beamwidth (symmetrical)
    %piston_angles=[0:0.1:20];
    %pist_model_best = beam_pattern(2*EPR(j),piston_angles,s2,fsint);
    %BW_3db(j) = 2*interp1(pist_model_best,piston_angles,-3);
    
    % Update waitbar
     waitbar(j/size(SIGNAL,1),h);
end
close(h)

figure(99), clf, box on, hold on
[B,BINT,RES,RINT,STATS] = regress(EPR(:),[range ones(length(range),1)]);
x = [5:0.5:20];
[P,S] = polyfit(range,EPR',1) ;
[y,delta] = polyconf(P,x,S,'predopt','curve'); %This is curve CI
color=0.95*[1 1 1];
h1 = fill([x fliplr(x)],[y-delta fliplr(y+delta)],color,'LineStyle','none','FaceColor',color,'FaceAlpha',1);
h2 = plot([x ; x]',[y-delta ; y+delta]','k--','color',[0.2 0.2 0.2]);
h3 = plot(x,y,'k','LineWidth',2) ;
h4 = plot(range,EPR,'ks','MarkerSize',3,'MarkerfaceColor',[0 0 0]);
ylabel('EPR (cm)','FontSize',14)
xlabel('Range (m)','FontSize',14)
set(gca,'XLim',[4 21],'YLim',[3 8],'XTick',[0:5:25])
disp(['RANGE R2: ' num2str(STATS(1)) ', F: ' num2str(STATS(2)) ', p = ' num2str(STATS(3))])
