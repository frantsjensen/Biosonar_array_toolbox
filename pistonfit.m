function [] = pistonfit(file,maxrange)
% function [] = pistonfit(file,maxrange)
% Estimate directional characteristics for on-axis biosonar clicks recorded
% on multi-element vertical hydrophone array using a parametric fit based
% on a circular piston
%
% Input:
% - file        file name (if file is in current directory) or 
%               path and file name of compiled click parameter file
%               (output from compileclickparams.m)
% - maxrange    Maximum range criterion for clicks (typically 20m)
%
% Output:
% - an output file with name given as outputfile = [inputfile '_beam.mat']
%   will be saved as the script is running. At first, this will contain
%   individual fit parameters (best estimate of EPR, and modelled beam
%   pattern) and after bootstrap procedure, this file will contain also
%   bootstrap information
%
% Function: 
% First, the composite beam pattern will be estimated through a single 
% parametric fit to the click amplitudes recorded across all hydrophones.
% Then, user can choose to continue estimating bootstrap confidence
% intervals around the estimated directional characteristics
%
% Bootstrapping is time-consuming and will lock up computer for hours. 
% An estimated time will be calculated. Make sure to plug in computer 
% while performing bootstrap.
%
% At the end of the bootstrap procedure, a fitted beam pattern figure
% and summary statistics (EPR, DI, -3 dB and -10 dB Beamwidth) will be
% provided. You can see summary statistics again by using showparams.m
%
% F. H. Jensen, 2013 (frants.jensen@gmail.com)

if nargin<2
    maxrange = 20 ;
end

if nargin<1
    file=[];
end

% If file is not defined, allow user to choose now
if isempty(file)
    [filename,pathname]=uigetfile('*.mat','Select compiled click parameter file');
    if filename==0, return, end % End script if user chooses cancel button
    file = fullfile(pathname,filename) ;
end

D=load(file);

if ~isfield(D,'ANGLES')
    error(' Selected file needs to include corrected angle of incidence')
end

% Find clicks measured within range criterion
k = find (D.RANGE < maxrange ) ;
ANGLES = D.ANGLES(k,:);
ASL = 10.^(D.RELASL(k,:)/10);   % Needs to be D.RELASL for new toolbox
N = size(ANGLES,1);             % N is number of on-axis clicks

% Find model on-axis click - signal with highest back-calculated amplitude
kmax    = min(find(D.ASLpp==max(D.ASLpp))) ;
s       = D.SIGNAL(kmax,:) ;
sr      = D.fs(kmax) ;

figure(99), plot([1:length(s)]/sr,s), title('supposed on-axis click')

userchoice = centermenu(' Accept default on-axis click? ','Accept','Choose wav file');
if userchoice == 2,      
    [filename,pathname]=uigetfile('*.wav','Select ideal on-axis wav file: ');
    if filename==0, return, end % End script if user chooses cancel button
    wavfile = fullfile(pathname,filename) ;
    [s,sr]=wavread(wavfile) ;
end

% Interpolate with factor 8
ip  = 8;                    % Interpolation factor
s   = interpolate2(s,ip);   % Interpolate signal
nsr =ip*sr;                 % New sample rate of interpolated signal

% Define piston diameters to test
pistsize=1:0.01:20;      % test these aperture diameters (cm)

% Do parametric fit for all piston sizes
angles=ANGLES(:); levels=ASL(:);

tic;
hwaitbar = waitbar(0,' Fitting piston to amplitude data ');
for  j=1:length(pistsize)
    pist = zeros(length(angles),1);
    for k=1:length(angles)
        pist(k)=energy(conv(imppist(abs(angles(k)),pistsize(j)/100,1524,nsr),s)); 
    end
    pist=pist/max(pist); %normalize
    
    % Calculate goodness of fit as the SSE between observed levels and expected levels
    %SSE(j)=sum((pist'-levels).^2);                         %squared error
    SSE(j) = sum((10*log10(pist)-10*log10(levels)).^2);    %log-squared error
    if ~rem(j,10) % Update waitbar (occasionally, to save computing time)
        waitbar(j/length(pistsize),hwaitbar)
    end
end
close(hwaitbar)

% Find time per iteration
timeperfit = toc; % We will use this to estimate time for bootstrap

% Now find piston size that minimizes error
err=sqrt(SSE)/sqrt(length(angles));         % normalize errors
[g merr]=min(err);                          % Find the minimum error

% EPR is half the diameter (cm)
EPR = 0.5*pistsize(merr) ;

%do this loop to get a nice plot at a larger angle range for the best piston size
piston_angles=[0:0.1:45];
for j=1:length(piston_angles) 
    pist_model(j)=energy(conv(imppist(piston_angles(j),pistsize(merr)/100,1524,nsr),s));
end
pist_model=pist_model/max(pist_model); %normalize to max

% The piston fit will be a little bit jagged, so we'll apply gentle smoothing
piston_model_best=smooth(piston_angles,10*log10(pist_model),0.1,'rloess');

% Create quick description
description = 'EPR in cm, angles in deg, piston model in dB energy rel peak energy';


% Make a couple default plots
figure(1), clf, % SSE as a function of piston size
subplot('position',[0.10 0.7 0.3 0.22]), hold on, box on
plot(pistsize,err,'k','LineWidth',2)
plot(pistsize(merr)*[1 1],[0 min(err)],'r','LineWidth',3)
xlabel('piston diameter (cm)','FontName','Helvetica','fontsize',12);
ylabel('SSE','FontName','Helvetica','fontsize',12);
    
% Best-fit piston compared to measured points
subplot('position',[[0.55 0.7 0.35 0.22]]), hold on, box on
plot(abs([-fliplr(piston_angles) piston_angles]),[fliplr(pist_model) pist_model],'k-',abs(-fliplr(angles)),fliplr(levels),'rx')
set(gca,'xlim',[0 30])
xlabel('Angle (degrees)','FontName','Helvetica','fontsize',12)
ylabel('Rel. energy','FontName','Helvetica','fontsize',12)

% Polar plot of best fitting piston model
subplot('position',[0.1 0.0 0.35 0.6]),
h=polar(-abs([-fliplr(piston_angles) piston_angles])/180*pi+pi/2,10*log10([fliplr(pist_model) pist_model])+40,'k');
set(gca,'ylim',[0 40],'xlim',[0 40]), set(h,'LineWidth',3)

% Piston fit and smoothed piston fit
subplot('position',[[0.55 0.1 0.35 0.4]]), hold on, box on, grid on
plot(piston_angles,10*log10(pist_model),'rx')
plot(piston_angles,piston_model_best,'k')
xlabel('Angle (degrees)','FontName','Helvetica','fontsize',12)
ylabel('Rel. energy','FontName','Helvetica','fontsize',12)

BW_3db = 2*interp1(piston_model_best,piston_angles,-3);
DI = 20*log10(185/BW_3db);
disp(['3 dB BW: ' num2str(BW_3db) 'degrees and DI: ' num2str(DI)])
clear BW_3db DI

% Now save data
dot = findstr(file,'.');
if dot
    outputfile = [file(1:dot-1) '_beam.mat'];
else
    outputfile = [file '_beam.mat'];
end
save(outputfile,'EPR','piston_angles','piston_model_best','ANGLES','ASL','N','description')
disp([ ' Saved beam fit (w/o bootstrap) to ' outputfile])



disp(' Now do bootstrap procedure to estimate confidence intervals around piston fit ')
disp(' If you have access to MatLab parallel computing toolbox, consider opening ')
disp(' worker pool using matlabpool.m and switching bootstrap for loop')
disp(' with a parfor loop')

userchoice = centermenu([' Bootstrap time estimated at ' num2str(2000*timeperfit/3600) ' hours'],' Do bootstrap now ',' Abort bootstrap ');

if userchoice == 1,

    % Now find the bootstrap distribution of the piston size
    B = 2000;
    angle_all  = ANGLES;    % We need to do bootstrap on a per-click basis
    levels_all = ASL ;      % to keep hydrophones together
    N = size(angle_all,1);  % N is number of on-axis clicks

    hwait = waitbar(0,'Bootstrapping');
    EPRboot = zeros(B,1);  % Save a bit of time by pre-allocating
    for nboot = 1:B, % For B bootstrap replicates
        % Select random clicks with replacement of indexes iboot
        iboot = ceil(N*rand(N,1)); 
        angle = angle_all(iboot,:); 
        levels = levels_all(iboot,:);
        angle=angle(:); levels=levels(:); 

        SSE = zeros(length(pistsize),1); % Preallocate
        for  j=1:length(pistsize)
            pist = zeros(length(angle),1) ;
            for k=1:length(angle)
                pist(k)=energy(conv(imppist(abs(angle(k)),pistsize(j)/100,1524,nsr),s)); 
            end
            pist=pist/max(pist); %normalize
            % Calculate goodness of fit as the SSE between observed levels and expected levels
            SSE(j) = sum((10*log10(pist)-10*log10(levels)).^2); %TESTING ONLY
        end

        % Find best fitting EPR
        err=sqrt(SSE)/sqrt(length(angle)); %normalize errors
        [g merr]=min(err); %wheres the minimum?
        EPRboot(nboot) = 0.5*pistsize(merr); % In cm
        % Update wait bar
        waitbar(nboot/B,hwait)
    end
    close(hwait)

    % Find 95% confidence interval: Two methods:

    % Only when normal distribution:
    % A.C. Davison and D.V. Hinkley (1996), p198-200
    alpha = 0.05;
    stat = EPR;
    bstat = EPRboot; % bootstrap statistics
    se = std(bstat);   % standard deviation estimate
    bias =mean(bstat-stat); % bias estimate
    za = norminv(alpha/2);   % normal confidence point
    lower = stat + bias + se*za; % lower bound
    upper = stat + bias - se*za;  % upper bound
    ci_norm = [lower;upper];  

    % percentile bootstrap CI
    pct1 = 100*alpha/2;
    pct2 = 100-pct1;
    lower = prctile(bstat,pct1); 
    upper = prctile(bstat,pct2);
    ci_prctile =[lower;upper];

    % Plot results and evaluate
    figure(2), clf, hold on, box on
    hist(EPRboot,20)
    hb(1)=plot([1 1]*stat,[0 0.5]*1000,'r','LineWidth',2);
    hb(2)=plot([1 1]*mean(bstat),[0 0.5]*1000,'r--','LineWidth',2);
    hb(3)=plot([1 1]*ci_prctile(1),[0 0.5]*1000,'g','LineWidth',2);
    plot([1 1]*ci_prctile(2),[0 0.5]*1000,'g','LineWidth',2)
    hb(4)=plot([1 1]*ci_norm(1),[0 0.5]*1000,'b--','LineWidth',2)
    plot([1 1]*ci_norm(2),[0 0.5]*1000,'b--','LineWidth',2)
    xlabel('Equivalent piston radius, cm','FontName','Helvetica','FontSize',12)
    ylabel('Bootstrap values','FontName','Helvetica','FontSize',12)
    
    legend(hb,'Best fit EPR','Mean bootstrap EPR','95% BCI normal','95% BCI percentile','Location','best')
    
    EPR_CI = ci_prctile;
    
    % Find click beam patterns
    piston_angles=[0:0.1:45];
    pist_model_best = beam_pattern(2*EPR,piston_angles,s,nsr);
    pist_model_lower = beam_pattern(2*EPR_CI(1),piston_angles,s,nsr); % Lower 95% CI
    pist_model_upper = beam_pattern(2*EPR_CI(2),piston_angles,s,nsr); % Upper 95% CI

    % Calculate -3 dB beamwidth (symmetrical)
    BW_3db = 2*interp1(pist_model_best,piston_angles,-3);
    BW_3db_CI = 2*[interp1(pist_model_upper,piston_angles,-3) interp1(pist_model_lower,piston_angles,-3)];

    % Calculate -10 dB beamwidth (symmetrical)
    BW_10db = 2*interp1(pist_model_best,piston_angles,-10);
    BW_10db_CI = 2*[interp1(pist_model_upper,piston_angles,-10) interp1(pist_model_lower,piston_angles,-10)];

    % Estimate Directivity Index (Following Madsen and Wahlberg 2007 approximation)
    DI = 20*log10(185/BW_3db) ;
    DI_CI = 20*log10(185./BW_3db_CI) ;
    
    % Save output of bootstrap procedure
    save(outputfile,'EPR','EPRboot','piston_angles','pist_model_best','pist_model_lower','pist_model_upper','BW_3db','BW_3db_CI','BW_10db','BW_10db_CI','DI','DI_CI','ANGLES','ASL','N','description')
    disp([ ' Saved beam fit with bootstrap to ' outputfile])

    fprintf('\n VARIABLE \t\t\t Mean \t [95%% BCI] \t\t units \n')
    fprintf(' EPR \t\t\t\t %2.2f \t [%2.2f-%2.2f] \t cm \n',EPR,EPR_CI(1),EPR_CI(2)) ;
    fprintf(' -3 dB beamwidth \t %2.2f \t [%2.2f-%2.2f] \t degrees \n',BW_3db,BW_3db_CI(1),BW_3db_CI(2)) ;
    fprintf(' -10 dB beamwidth \t %2.2f \t [%2.2f-%2.2f] \t degrees \n',BW_10db,BW_10db_CI(1),BW_10db_CI(2)) ;
    fprintf(' DI \t\t\t\t %2.1f \t [%2.1f-%2.1f] \t dB \n',DI,DI_CI(2),DI_CI(1)) ;
    fprintf(' N  \t\t\t\t %2.0f \t \n\n',N) ;

    figure(3), clf, hold on, box on,
    h2 = plot(abs(ANGLES),10*log10(ASL),'ks','MarkerFaceColor','k','MarkerSize',3);
    h2(2) = plot(piston_angles,pist_model_best,'k','LineWidth',3,'color',[0.3 0.3 0.3]);
    h2(3) = plot(piston_angles,pist_model_lower,'k--','LineWidth',2,'color',[0.3 0.3 0.3]);
    h2(4) = plot(piston_angles,pist_model_upper,'k--','LineWidth',2,'color',[0.3 0.3 0.3]);
    set(gca,'fontsize',12,'ylim',[-25 1],'xlim',[0 25])
    xlabel('Angle of incidence (degrees vertical)','FontName','Helvetica','FontSize',14)
    ylabel('Power (dB re. max)','FontName','Helvetica','FontSize',14)
    ht=text(14,-3,sprintf('-3 dB beamwidth \n Mean: %2.1f \n 95%% BCI: %2.1f to %2.1f ',[BW_3db,BW_3db_CI(1),BW_3db_CI(2)]),'FontName','Helvetica','FontSize',12);
    legend(h2(1:3),'Array data','Piston fit','Piston fit bootstrap CI','Location','SouthWest')
    legend boxoff
end