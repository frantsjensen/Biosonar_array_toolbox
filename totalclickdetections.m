function [numclicks] = totalclickdetections(filepath)
% [numclicks] = totalclickdetections(filepath)
% Looks through all click detection files in folder
% specified in filepath to find ICI and total number of
% clicks detected. Uncorrected ICI statistics are for all clicks detected.
% Corrected ICI are only for clicks that have been through
% visualizeclicktrain and where overlapping click trains will have been
% separated.
%
% If no filepath is given, current click detector folder is searched
% (i.e. subfolder /clickdetector/ in data directory specified in toolbox.m)
%
% F. H. Jensen, 2013 (frants.jensen@gmail.com)

SETTINGS = toolbox_settings();

if nargin<1,
    filepath = fullfile(SETTINGS.datadir,'clickdetector') ;
end

% Find data directory and list of files within it
D = dir (filepath);

% Preallocate
ICI = [];
ALLICI = [];
numclicks = 0;
n=0;

% Load in all data files one by one and accumulate list of parameters
for j = 1:length(D),
    if strfind ( D (j).name , 'clicks.mat') % Check that the file is a click detection file
        
        % Load in click parameters
        C = load (fullfile(filepath,getfield( D(j) , 'name')));
        
        n = n+1;
     
        % Source information
        numclicks(n,:)  = length(C.C);      % Number of detected clicks
        ici = C.ICI(find(ICI>0)) ;
        ALLICI = diff(C.C)/C.fs ;           % Corrected ICI (no overlapping click series)
        ICI = [ICI ; ici(:)] ;              % Uncorrected ICI
    end
end

numclicks = sum(numclicks) ;

% Write out summary to command prompt
fprintf('\n Total click detections \t\t %5.0f \n',numclicks) ;
fprintf(' Corrected ICI \t\t\t\t %3.1f +/- %3.1f \t\t [%3.1f-%4.1f] \t ms \n',easystats(1000*ICI)) ;
fprintf(' Uncorrected ICI \t\t\t %3.1f +/- %3.1f \t\t [%3.1f-%4.1f] \t ms \n',easystats(1000*ALLICI)) ;
