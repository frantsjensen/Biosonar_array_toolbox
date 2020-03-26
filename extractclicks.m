% EXTRACTCLICKS - find and analyze on-axis biosonar clicks
%
% extractclicks (file) works through N-channel wave file with path and 
% file name specified in 'file' to detect, localize and plot biosonar 
% clicks, and allow user to interactively select a subset of biosonar 
% clicks for quick estimation of source parameters.
%
% extractclicks (file,gain) works as above but replaces the amplification
% gain specified in toolbox_settings.m with the gain specified as input,
% and then recalculates clip levels for each hydrophone.
%
% extractclicks (file,gain,threshold) works as above but uses a
% user-defined (constant) threshold to detect clicks. Use gain = [] to
% preserve gain defined in toolbox_settings file.
%
% TOOLBOX_SETTINGS.M: All properties of this and the underlying scripts
% are regulated by toolbox_settings.m - please adjust this before running
% toolbox. A copy of current toolbox settings is saved with all data.
% 
% OUTPUT: Click detections and individual click source parameters are 
% saved to user defined folder (click parameters) and subfolder 
% (click detections). 
%
% Once all data has been analysed and all on-axis clicks have been 
% extracted, run compileclickparams.m to accumulate all click parameters 
% into one file, then run pistonfit.m to estimate vertical beam pattern 
% and bootstrap confidence intervals.
%
% F. H. Jensen, 2013 (frants.jensen@gmail.com)

%%%%%%%%%%%%%%% A: Define Parameters and Files %%%%%%%%%%%%%%%
function [] = extractclicks (file,gain,threshold)

% Insert dynamic threshold
if nargin < 3,
    threshold=[];
end

% Insert empty placeholder for signal gain
if nargin <2,
    gain = [];
end

% Check file name
if ~ischar(file)
    error('Need a text string as input - see help')
elseif ~exist(file,'file')
    error('A file with the given name and path does not exist')
elseif ~strcmp(file(end-3:end),'.wav')
    error('Input file needs to be a wave file ending with .wav')
end

% Replace '\' with '/' (mac compatibility)
file(strfind(file,'\'))='/';

% Find file path and file name
fpath = file(1:max(strfind(file,'/')));
fname = file(max(strfind(file,'/'))+1:end);

% Adjust toolbox settings
SETTINGS = toolbox_settings();

% Now adjust for different gain compared to settings file
if ~isempty(gain)
    SETTINGS.clip_level = SETTINGS.clip_level +SETTINGS.gain - gain ;
    SETTINGS.gain = gain ;
    disp([' Gain adjusted to ' num2str(gain) ' dB compared to settings file']);
end

%%%%%%%%%%%%%%% A: Detect Clicks %%%%%%%%%%%%%%%

% ALTERNATIVE A1: Manually adjust threshold based on peaks of envelope (clicks)
% [time,maxenv]=findthreshold(file,SETTINGS);

% %Plot maximum envelope (peak detection)
% figure(1), clf, plot(time,maxenv,'k')
% set(gca,'YScale','log')
% threshold = input('Enter Detection Threshold (constant detection level best): ');
% [C,fs,newthreshold] = detectclicks (file , threshold , SETTINGS);

% ALTERNATIVE A2: CLICK DETECTION THRESHOLD BASED ON NOISE LEVEL OF RECORDING
[C,fs,newthreshold] = detectclicks (file ,threshold, SETTINGS);
if length(newthreshold)>1,
    disp(['Threshold 5th, 50th, 95th percentile: ' num2str(prctile(newthreshold,[5 50 95])) ])
end

%%%%%%%%%%%%%%% B: Locate Clicks %%%%%%%%%%%%%%%
Sxy   = zeros(length(C),2);
range = zeros (length(C),1);
RLpp  = zeros (length(C),length(SETTINGS.R));
sigcenter = zeros (length(C),SETTINGS.NFFT) ;

% Go through clicks and localize each
hwaitbar = waitbar(0,' Localizing detected clicks ');
for click_j = 1:length(C),
    [LOC,RLtemp,Tc,SIGTEMP] = localizeclick (file,C(click_j),SETTINGS);
    Sxy(click_j,:) = LOC;
    range(click_j)  = sqrt(LOC(1)^2+(LOC(2)-SETTINGS.R(SETTINGS.CH))^2); % Dist to hydrophone used for click detection
    RLpp (click_j,:) = RLtemp;
    sigcenter (click_j,:) = SIGTEMP(:,SETTINGS.CH)' ;
    waitbar(click_j/length(C),hwaitbar)
end
close(hwaitbar)

% Create quick figure of all click detections
figure(2), clf,
subplot(3,1,1), plot(C/fs,real(Sxy(:,1)),'ro'); ylabel('Distance, m'), xlabel('Time, s')
subplot(3,1,2), plot(C/fs,max(RLpp,2),'ko'); ylabel('RL, dB'), xlabel('Time, s')
subplot(3,1,3), plot(C(2:end)/fs,diff(C)/fs,'ko'); ylabel('ICI, s'), xlabel('Time, s'), set(gca,'YLim',[0 ceil(5*median(diff(C)/fs))])

% Save preliminary click detector output
if ~exist(fullfile(SETTINGS.datadir,'clickdetector'))
    mkdir(SETTINGS.datadir,'clickdetector')
    disp([' New folder created: ' fullfile(SETTINGS.datadir,'clickdetector')])
end
savefile = fullfile(SETTINGS.datadir,'clickdetector',[fname(1:end-4) '_clicks.mat']);
save (savefile, 'file','fs','C','RLpp','newthreshold','range','SETTINGS','-v6');

% Output user message
fprintf('\n Data from click detector saved \n\n')

%%%%%%%%%%%%%%% C: Identify On-axis Clicks %%%%%%%%%%%%%%%
[ICI,Cx] = visualizeclicktrain(C,RLpp,range,file,SETTINGS);

% Save click detector output with on-axis click information
save (savefile, 'file','fs','C','RLpp','newthreshold','range','ICI','Cx','SETTINGS','-v6');

% Output user message
fprintf('\n Data from click detector saved with information on on-axis click location \n\n')