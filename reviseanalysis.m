function [] = reviseanalysis(filepath)

% function [] = reviseanalysis(filepath)
%
% Finds the source parameters of biosonar clicks that have already been 
% extracted from sound files, using information about source files and 
% click locations contained in individual *params.mat files.
%
% reviseanalysis(filepath) will modify all click parameter files located in
% the specified folder. Reviseanalysis() will modify all click parameter
% files located in the default data directory specified in the current
% toolbox settings.
%
% reviseanalysis will call the modified click analysis script 
% analyzeclick_revision.m for each existing *params.mat file, and save
% all newly calculated click parameters directly to the same parameter file.
% Thus, these two scripts will allow you to change analysis parameters 
% (such as FFT size) or add analysis components (such as different metrics)
% after clicks have already been extracted.
%
% Change the subscript, analyzeclick_revision.m, to reflect new analysis
% procedure.
%
% All parameters extracted following equations in Madsen and Wahlberg 2007
% and W. W. L. Au (1993): Sonar of Dolphins
%
% F. H. Jensen, 2014 (frants.jensen@gmail.com)

% Set filepath to default filepath defined in toolbox_settings.m
SETTINGS = toolbox_settings();

if nargin<1,
    filepath = SETTINGS.datadir ;
end

% Find data directory and list of files within it
D = dir (filepath);
n=0;

% Load in all data files one by one and accumulate list of parameters
for i = 1:length(D),
    if strfind ( D (i).name , 'params.mat') % Check that the file is a click parameter file
        
        % Make settings structure that can be updated with clip level
        NEWSETTINGS=SETTINGS;
        
        % Find path and name of current click parameter file
        currentfile = [filepath getfield( D(i) , 'name')] ;
        
        % Load in click parameters
        C = load (currentfile);
        
        % We need some of the original data
        file            = C.file;   % Original source file
        Tc              = C.Tc;     % peak sample on each receiver
        LOC             = C.SOURCE_LOC; % source location
        ICI             = C.ICI;
        
        % Check if settings have different clip level than when the click
        % was analysed
        if ~all(SETTINGS.clip_level==C.SETTINGS.clip_level)
            disp([' File: ' getfield(D(i),'name') ' Clip level difference: '])
            disp([' Current clip level (' num2str(SETTINGS.clip_level(1)) ' dB) vs previous clip level ( ' num2str(C.SETTINGS.clip_level(1)) ' dB)'])
            disp([' Using previously configured clip level '])
            NEWSETTINGS.clip_level=C.SETTINGS.clip_level;
            NEWSETTINGS.gain=C.SETTINGS.gain;
        end
        
        % Then we pass that original data to the revised analyze click script
        % Changes will be saved directly to old click parameter file
        analyzeclick_revision(file,Tc,LOC,ICI,NEWSETTINGS,currentfile) ;
        n=n+1;
    end
end

disp([ ' A total of ' num2str(n) ' click parameter files in ' filepath ' updated'])