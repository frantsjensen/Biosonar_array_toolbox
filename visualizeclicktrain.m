function [ICIcorr,Cx] = visualizeclicktrain(C,RL,range,file,SETTINGS)
% VISUALIZECLICKTRAIN(C,RLpp,range,fs,file,SETTINGS) separates overlapping
% click series from each other and plots signal received level and 
% estimated range for each click in a scan, after which it allows the user
% to select on-axis clicks for further source parameter estimation
%
% Input:
%   C:  	click peak sample for each click
%   RL:     received levels across all hydrophones for each click
%   range:  horizontal range for each click
%   file:   string containing full path and file name to source file
%   SETTINGS is a structure from toolbox_settings.m
%
% Selection of on-axis clicks:
% A plot of detected clicks for each scan will be shown. Left-clicking
% within that plot marks a click and allows user to estimate and save
% source parameters for that click. To proceed to next scan, press 'f'
% (forward).
%
% Interactive buttons:
% Press 'f' or 'F' to go to next scan
% press 't' to zoom in (x2) near cursor
% press 'y' to zoom out and center on cursor
% press 'r' to reset to original zoom
% press 's' to save all ICIs for red clicks (corrected ICI stored in 
%              click detector files *clicks.mat)
% left click over on-axis biosonar click to extract parameters
%
%
% [ICI,Cx] = visualizeclicktrain(C,RLpp,range,fs,file,SETTINGS)
% will return the corrected ICI for each click in a scan (or 0 otherwise)
% and the sample corresponding to each on-axis click selected by user
%
% F. H. Jensen, 2013 (frants.jensen@gmail.com)

% Load settings
maxici      = SETTINGS.max_ICI;
minclicks   = SETTINGS.min_clicks ;
R           = SETTINGS.R ;
symbols     = SETTINGS.symb ;
siz         = SETTINGS.siz ;

% Create temporary results
ICI = zeros(1,length(C));   % ICI for each click w click train separation
ICIcorr = zeros(1,length(C));

% Check if any clicks have already been analyzed
[Cx,Cx_file,Cx_RLpp] = findanalyzedclicks(file,SETTINGS);

% Find sample rate
[temp , fs] = wavexread ( file , [1 2] );

% Divide the entire click sequence into scans based on maximum ICI
ICIt = 1000*[0 diff(C)]/fs; 
scan_index = [1 find(ICIt>maxici) length(ICIt)];

% Find scans longer than minimum clicks setting
acceptedscans = find(diff(scan_index)>minclicks);

% Start at first acceptable scan
scan_k = min(acceptedscans);

% Look through scans
while ~isempty(scan_k) ;

    % Find clicks within this scan
    k=[scan_index(scan_k):scan_index(scan_k+1)-1];

    % Find axes for plot
    axy = [floor(min(min(RL(k,:)))) ceil(max(max(RL(k,:))))];

    % Separate overlapping click trains based on ICI
    groups = separateclickseries(C(k)) ;
    k1 = k(find(groups==1));
    k2 = k(find(groups==2));

    % Recalculate corrected ICI for main click series)
	%ICIcorr(k1) = [0 diff(C(k1))]'/fs ; %BUG! THIS LINE ACTIVE IN CURRENT
	%DISTRIBUTED CODE
            
    % For each hydrophone, plot clicks in this scan in chosen color pattern
    figure(10), clf, set(gcf,'Name','Scan visualization')
    hax(1)=axes('position',[0.12 0.4 0.8 0.5]); hold on, box on
    for j=1:length(R), 
        h1(j) = plot(C(k1)/fs,RL(k1,j),['k' char(symbols(j))],'MarkerSize',siz(j),'markerfacecolor',[1 0 0]*sin(pi*j/(length(R)+1)),'markeredgecolor',[1 0 0]*sin(pi*j/(length(R)+1)) );
        legg{j} = ['Ch. ' num2str(j)] ;
        if ~isempty(k2)
            h2(j) = plot(C(k2)/fs,RL(k2,j),['k' char(symbols(j))],'MarkerSize',siz(j),'markerfacecolor',[0 0 1]*sin(pi*j/(length(R)+1)),'markeredgecolor',[0 0 1]*sin(pi*j/(length(R)+1)));
        end
    end
    tit=['Scan No. ' num2str(scan_k) ' of ' num2str(length(scan_index)-1) ' showing clicks ' num2str(k(1)) ' to ' num2str(k(end))];
    title(tit,'FontName','Helvetica','FontSize',14,'FontWeight','Bold', 'HorizontalAlignment' , 'Center')
    xlabel('Time into file, s','FontName','Helvetica','FontSize',14)
    ylabel('Received level, dB','FontName','Helvetica','FontSize',14)
    axis([0.1*floor(10*C(k(1))/fs) 0.1*ceil(10*C(k(end))/fs) axy(1) axy(2)+3]);
    hleg =legend(h1,legg,'location','best'); set(hleg,'box','off')

    hax(2)=axes('position',[0.12 0.1 0.8 0.2]); hold on, box on
    if ~isempty(k2)
        plot(C(k2)/fs,range(k2),'ks','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0.5]);
    end
    plot(C(k1)/fs,range(k1),'ko','MarkerEdgeColor',[0.5 0 0],'MarkerFaceColor',[1 0 0]);
    if min(range(k))<100,
        axis([0.1*floor(10*C(k(1))/fs) 0.1*ceil(10*C(k(end))/fs) max([0 floor(min(range(k)))]) min([100 3+ceil(max(range(k)))]) ])
    else
        axis([ 0.1*floor(10*C(k(1))/fs) 0.1*ceil(10*C(k(end))/fs) max([0 floor(min(range(k)))]) 3+ceil(max(range(k))) ])
    end
    xlabel('Time into file, s','FontName','Helvetica','FontSize',14)
    ylabel('Range, m','FontName','Helvetica','FontSize',14)

	% Plot all analyzed on-axis clicks
    axes(hax(1)) ;
    xlim = get(gca,'XLim') ;
    for Cx_k = 1:length(Cx)
        if all( [ xlim(1) < Cx(Cx_k)/fs ; Cx(Cx_k)/fs < xlim(2)] )
            kclose = [max(find(Cx(Cx_k)>C)) min(find(Cx(Cx_k)<C))] ;
            dist = abs(C(kclose)/fs-Cx(Cx_k)/fs) ;
            kclose = kclose(find(dist==min(dist)));
            xm = Cx(Cx_k)/fs + ICIt(kclose)*[-0.5 0.5]/1000 ;
            ym = Cx_RLpp(Cx_k) ;
            plot([xm(1) xm(1) NaN xm(1) xm(2) NaN xm(2) xm(2)],[0.5 1.5 NaN 1 1 NaN 0.5 1.5]+ym,'k','LineWidth',3);
            text(mean(xm),ym+2,Cx_file(Cx_k),'verticalalignment','bottom','horizontalalignment','center')
        end
    end
    
    % Link the x-axis so that zooming in on one plot also zooms in on
    % second plot, at least in x-axis
    linkaxes(hax,'x');
    
    done = 0 ;
    while done == 0,
        axes(hax(1)) ; htemp=plot(1,1,'kx'); delete(htemp)
                
        [gx gy button] = ginput(1) ;
        
        if button=='f' | button=='F'
            done = 1 ;
            
        elseif button=='t' % Zoom in factor 2
            lim     = axis;
            outside = gx<lim(1) || gx>lim(2) || gy<lim(3) || gy>lim(4);
            if ~outside,
            	% Center on click (x-axis only) and expand x-axis by x 0.5
                lim = [gx+diff(lim(1:2))/2*[-1 1]*0.5 lim(3:4)];
                axis(lim)            
            end
                
            
        elseif button =='y' % Zoom out factor 2
            lim     = axis;
            outside = gx<lim(1) || gx>lim(2) || gy<lim(3) || gy>lim(4);
            if ~outside,
            	% Center on click (x-axis only) and expand x-axis by x 2
                lim = [gx+diff(lim(1:2))/2*[-1 1]*2 lim(3:4)];
                axis(lim)            
            end
            
        elseif button =='r';
            set(gca,'xlim',[0.1*floor(10*C(k(1))/fs) 0.1*ceil(10*C(k(end))/fs)]);
                
            
        elseif button =='s' % Manually re-categorize click to different series
             % Recalculate corrected ICI for main click series)
            ICIcorr(k1) = [0 diff(C(k1))]'/fs ;
            disp([num2str(length(k1)) 'Corrected ICIs saved'])
            
        elseif button==1 | button=='p'
            kclose = [max(find(gx>C/fs)) min(find(gx<C/fs))] ;
            dist = abs(C(kclose)/fs-gx) ;
            kclose = kclose(find(dist==min(dist)));
            
            % Create a marker for this click
            xm = C(kclose)/fs + ICIt(kclose)*[-0.5 0.5]/1000 ;
            ym = ceil(max(RL(kclose,:))) ;
            hm = plot([xm(1) xm(1) NaN xm(1) xm(2) NaN xm(2) xm(2)],[0.5 1.5 NaN 1 1 NaN 0.5 1.5]+ym,'k','LineWidth',3);
            
            userchoice = centermenu(' Choose what to do with this click ','Extract parameters','Abort');
            
            if userchoice==1
                % Do source localization, then extract and save parameters
                [LOC,RLtemp,Tc] = localizeclick(file,C(kclose),SETTINGS,1) ;
                [clickname]=analyzeclick(file,Tc,LOC,ICIt(kclose),SETTINGS) ;
                
                if ~isempty(clickname)
                    % Add text to original plot
                    axes(hax(1)); text(mean(xm),ym+2,clickname,'verticalalignment','bottom','horizontalalignment','center')

                    % Remember location of on-axis click for host program
                    Cx = [Cx ; C(kclose)];
                    Cx_file{length(Cx_file)+1} = clickname ;
                    Cx_RLpp = [Cx_RLpp ; max(RLtemp)] ;
                else
                    delete(hm)
                end
            else
                delete(hm)
            end
            
            
            
        % QUIT CLICK VISUALIZATION PROGRAM
        elseif button=='q',
            disp(['Exiting click visualization tool'])
            return
            
        end
    end
    
    if scan_k == max(acceptedscans),
        userchoice = centermenu(' End of file ','Exit click train visualization','Return to beginning of file');
        if userchoice == 2,
            scan_k = 0 ; % Reset to 0 and find first acceptable scan three lines down
        end
    end
    scan_k = acceptedscans(min(find(acceptedscans>scan_k)));
end
    
%%%%%% subfunction to find already analyzed clicks %%%%%%

function [Cx,Cx_file,Cx_RLpp] = findanalyzedclicks(file,SETTINGS)

filepath = SETTINGS.datadir ;
ch = SETTINGS.CH ;
Cx = []; Cx_file = {} ; Cx_RLpp = [];

D = dir(filepath) ;
for k=1:length(D),
    if strfind ( D (k).name , 'params.mat') % Check that the file is a click parameter file
        
        % Load in click parameters
        C = load ([filepath getfield( D(k) , 'name')]);
        
        dirseps = sort([strfind(file,'/') strfind(file,'\')]) ;
        
        % Build up column vector of analyzed clicks
        if strfind(C.file,file(dirseps(end)+1:end))
            seps = strfind(D(k).name,'_') ;
            Cx = [Cx ; C.Tc(ch)] ; 
            Cx_RLpp = [Cx_RLpp ; max(C.RLpp)] ;
            Cx_file{length(Cx_file)+1} = D(k).name(seps(end-1)+1:seps(end)-1) ;
        end
    end 
end