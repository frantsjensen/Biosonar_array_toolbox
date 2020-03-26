function [] = showclickparams()
% SHOWCLICKPARAMS allows user to interactively choose .mat files with 
% click parameter data and display summary statistics.
% User will be prompted to locate file with compiled click parameters
% (output from compileclickparams.m) and summary statistics will be shown.
% Then, user will be prompted to locate file with beam pattern parameters
% (from pistonfit.m including bootstrap data) and directivity information
% will be shown.
%
% F. H. Jensen, 2013 (frants.jensen@gmail.com)

% Find compiled click parameter file
[filename,pathname]=uigetfile('*.mat','Select compiled click parameter file: ');
if filename==0, return, end % End script if user chooses cancel button
file = fullfile(pathname,filename) ;
load(file);

% Show click parameter data
fprintf(' VARIABLE \t\t\t Mean+/-S.D. \t\t [min-max] \t\t units \n')
fprintf(' SLpp\t\t\t\t %3.1f +/- %2.1f \t\t [%3.1f-%3.1f] \t dB re 1 muPa pp @1m \n',easystats(ASLpp)) ;
fprintf(' SLrms (-10 dB)\t\t %3.1f +/- %2.1f \t\t [%3.1f-%3.1f] \t dB re 1 muPa rms  @1m \n',easystats(ASLrms_amp)) ;
fprintf(' SLefd (-10 dB)\t\t %3.1f +/- %2.1f \t\t [%3.1f-%3.1f] \t dB re 1 muPa rms  @1m \n',easystats(ASLflux_amp)) ;
fprintf(' D (-10 dB)\t\t\t %2.1f +/- %1.1f \t\t [%2.1f-%2.1f] \t mus \n',easystats(DUR_amp)) ;
fprintf(' Fc \t\t\t\t %3.1f +/- %2.1f \t\t [%3.1f-%3.1f] \t kHz \n',easystats(Fc)) ;
fprintf(' Fp \t\t\t\t %3.1f +/- %2.1f \t [%3.1f-%3.1f] \t kHz \n',easystats(Fp)) ;
fprintf(' BW (-3 dB) \t\t %3.1f +/- %2.1f \t\t [%3.1f-%3.1f] \t kHz \n',easystats(BW_3db)) ;
fprintf(' BW (-10 dB) \t\t %3.1f +/- %2.1f \t\t [%3.1f-%3.1f] \t kHz \n',easystats(BW_10db)) ;
fprintf(' BW (rms) \t\t\t %3.1f +/- %2.1f \t\t [%3.1f-%3.1f] \t kHz \n',easystats(BW)) ;
fprintf(' Q (rms)  \t\t\t %2.1f +/- %1.1f \t\t [%2.1f-%2.1f] \t\t - \n',easystats(Qrms)) ;
fprintf(' ICI \t\t\t\t %3.1f +/- %2.1f \t\t [%3.1f-%3.1f] \t ms \n',easystats(1000*ICI)) ;
fprintf(' RANGE \t\t\t\t %3.1f +/- %2.1f \t\t [%3.1f-%3.1f] \t m \n',easystats(RANGE)) ;
fprintf(' N \t\t\t\t\t\t %3.0f \n',length(ASLpp)) ;



% Now find beam pattern data
[filename,pathname]=uigetfile('*.mat','Select file with beam pattern data: ');
if filename==0, return, end % End script if user chooses cancel button
file = fullfile(pathname,filename) ;
load(file)

% Show beam pattern data
fprintf('\n VARIABLE \t\t\t Mean \t [95%% BCI] \t\t units \n')
fprintf(' EPR \t\t\t\t %2.2f \t [%2.2f-%2.2f] \t cm \n',EPR,EPR_CI(1),EPR_CI(2)) ;
fprintf(' -3 dB beamwidth \t %2.2f \t [%2.2f-%2.2f] \t degrees \n',BW_3db,BW_3db_CI(1),BW_3db_CI(2)) ;
fprintf(' -10 dB beamwidth \t %2.2f \t [%2.2f-%2.2f] \t degrees \n',BW_10db,BW_10db_CI(1),BW_10db_CI(2)) ;
fprintf(' DI \t\t\t\t %2.1f \t [%2.1f-%2.1f] \t dB \n',DI,DI_CI(2),DI_CI(1)) ;
fprintf(' N  \t\t\t\t %2.0f \t \n\n',N) ;


