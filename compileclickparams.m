function [datafile,sig,fs] = compileclickparams(filepath,dataname,maxrange)
% [datafile,sig,fs] = compileclickparams(filepath,dataname,maxrange)
% Compiles all click source parameters from parameter files (*params.mat)
% in the folder defined by the string 'filepath', filters according to 
% maximum range defined in 'maxrange', and saves to file name defined in 
% the string 'dataname'.
%
% Output:
% - datafile    Name of output .mat file containing parameters
% - sig         Waveform of signal with highest ASL
% - fs          Sample rate of signal with highest ASL
%
% Use wavwrite('speciesname.wav',sig,fs) to save model waveform 
% Use pistonfit.m to estimate directional characteristics
%
% F. H. Jensen, 2013 (frants.jensen@gmail.com)



% Set filepath to default filepath defined in toolbox_settings.m
SETTINGS = toolbox_settings();

% Set default range criterion
if nargin<2
    maxrange = 50 ; % Range criterium
    disp(['Compiling data with range less than ' num2str(maxrange)])
end

% Construct default output name
if nargin<3  % Output data file name, format: 'Author_species_year'
    dataname = [ SETTINGS.clickprefix '_' date]; % Default only prefix+date
end

if nargin<1,
    filepath = SETTINGS.datadir ;
end

% Find data directory and list of files within it
D = dir (filepath);

% Make counter of on-axis signals
n = 0;

% Load in all data files one by one and accumulate list of parameters
for i = 1:length(D),
    if strfind ( D (i).name , 'params.mat') % Check that the file is a click parameter file
        
        % Load in click parameters
        C = load ([filepath getfield( D(i) , 'name')]);
        
        % Skip this click if range is outside of confident localization range
        if C.RANGE(1)>maxrange
            continue
        end
        
        % Find channel with on-axis click
        [temp,ch]           = max(C.ASLrms_amp);        

        % Ignore clicks that are off-axis after compensating for range
        % (only important for very close range clicks)
        if ch==1 | ch>(length(C.ANGLE)-1)
            disp(D(i).name)
            continue
        end
                
        % Save data
        n = n+1;
        FOCUS(n,:)          = ch;
        
        paramsfile{n}             = D(i).name;
        
        % Source information
        file{n}             = C.file;   % Original source file
        sample(n,:)         = C.Tc;     % peak sample on each receiver
        cliplevel(n,:)      = C.SETTINGS.clip_level;
        fs (n,:)            = C.fs;
        
        % Repetition rate
        ICI(n,:)            = C.ICI;
        
        % Peak-peak amplitude
        ASLpp(n,:)          = C.ASLpp(ch);

        % -10 dB amplitude measure
        ASLrms_amp(n,:)     = C.ASLrms_amp(ch);
        ASLflux_amp(n,:)    = C.ASLflux_amp(ch);
        DUR_amp(n,:)        = C.WINDOW_amp(ch);

        % 95% energy measure
        ASLrms_en(n,:)      = C.ASLrms_en(ch);
        ASLflux_en(n,:)     = C.ASLflux_en(ch);
        DUR_en(n,:)         = C.WINDOW_en(ch);
        
        % Estimated range and source location
        RANGE(n,:)          = C.RANGE(ch)';
        Sxy(n,:)            = C.SOURCE_LOC(:)';

        % Spectral parameters
        Fp(n,:)             = C.Fp(ch);
        Fc(n,:)             = C.Fc(ch);       
        BW(n,:)             = C.RMSBW(ch);
        BW_3db(n,:)         = C.BW_3db(ch);
        BW_10db(n,:)        = C.BW_10db(ch);
        Qrms(n,:)           = C.Fc(ch)/C.RMSBW(ch);         
        Q_3dB(n,:)          = C.Fp(ch)./C.BW_3db(ch);

        % Save also all rms amplitudes for use in estimating DI:
        ASLall(n,:)         = C.ASLrms_amp(:)';     % dB re 1 muPa
        THETA(n,:)          = C.ANGLE(:)'; % degrees rel to channel ch
        
        % Estimate exact angle of incidence and relative ASL for DI estimate
        if size(C.SETTINGS.R,2)==1
            chs                 = min([max([ch 2]) length(C.ASLrms_amp)-1])+[-1 0 1];
            [thisangle thisamp] = lagrange3pt(THETA(n,chs),10.^((ASLall(n,chs)-120)/20));
            ANGLES(n,:)         = THETA(n,:)-thisangle;    
            RELASL(n,:)         = ASLall(n,:)-120-20*log10(thisamp);
            if any(RELASL(n,:)>0),
                disp('Warning: Relative ASL should not exceed 0 dB (perfectly onaxis)')
                disp(['File: ' D(i).name])
            end
            clear thisangle thisamp chs
        elseif size(C.SETTINGS.R,2)==2,
            error('2D arrays not yet supported')
            [angles,relasl]     = gridfit(C.ANGLE(:)',C.ASLrms_amp(:)',C.SETTINGS.R) ;
            ANGLES(n,:)         = angles ;
            RELASL(n,:)         = relasl ;
            clear angles relasl
        end
                
        % Save extracted synchronized waveforms and power spectral density
        SIGS{n}             = C.SIGNAL ;        % Raw waveform
        SPECTRA{n}          = C.SPEC ;          % Power spectral density
        
        % Save on-axis waveform and power spectral density
        SIGNAL(n,:)         = C.SIGNAL(:,ch)' ; % On-axis waveform
        if n>1, 
            if length(C.SPEC(:,ch))~=size(SPEC,2)
                disp(['Error in ' getfield( D(i) , 'name') ': Different fft size or interp'])
                continue
            end
        end
        SPEC(n,:)           = C.SPEC(:,ch)' ;   % On-axis PSD
        FREQ(n,:)           = C.FREQ ;          % Frequency axis for PSD
        
        clear C temp ch
    end
end

if n==0;
    error(['No click parameter files found in ' filepath])
end

% Clear values that should not be saved
clear D i SETTINGS filepath

% Save data to summary file
if ~isempty(findstr(dataname,'.'))
    dataname=dataname(1:findstr(dataname,'.')-1);
end
save ([dataname '.mat'])
datafile = [dataname '.mat'] ;
fprintf(['\n ' num2str(n) ' clicks processed and saved as ' datafile '\n\n'])

% Write out summary to command prompt
fprintf(' VARIABLE \t\t\t Mean+/-S.D. \t\t [min-max] \t\t units \n')
fprintf(' SLpp\t\t\t\t %3.1f +/- %2.1f \t\t [%3.1f-%3.1f] \t dB re 1 muPa pp @1m \n',easystats(ASLpp)) ;
fprintf(' SLrms (-10 dB)\t\t %3.1f +/- %2.1f \t\t [%3.1f-%3.1f] \t dB re 1 muPa rms  @1m \n',easystats(ASLrms_amp)) ;
fprintf(' SLefd (-10 dB)\t\t %3.1f +/- %2.1f \t\t [%3.1f-%3.1f] \t dB re 1 muPa rms  @1m \n',easystats(ASLflux_amp)) ;
fprintf(' D (-10 dB)\t\t\t %2.1f +/- %1.1f \t\t [%2.1f-%2.1f] \t mus \n',easystats(DUR_amp)) ;
fprintf(' Fc \t\t\t\t %3.1f +/- %2.1f \t\t [%3.1f-%3.1f] \t kHz \n',easystats(Fc)) ;
fprintf(' Fp \t\t\t\t %3.1f +/- %2.1f \t\t [%3.1f-%3.1f] \t kHz \n',easystats(Fp)) ;
fprintf(' BW (-3 dB) \t\t %3.1f +/- %2.1f \t\t [%3.1f-%3.1f] \t kHz \n',easystats(BW_3db)) ;
fprintf(' BW (-10 dB) \t\t %3.1f +/- %2.1f \t\t [%3.1f-%3.1f] \t kHz \n',easystats(BW_10db)) ;
fprintf(' BW (rms) \t\t\t %3.1f +/- %2.1f \t\t [%3.1f-%3.1f] \t kHz \n',easystats(BW)) ;
fprintf(' Q (rms)  \t\t\t %2.1f +/- %1.1f \t\t [%2.1f-%2.1f] \t\t - \n',easystats(Qrms)) ;
fprintf(' ICI \t\t\t\t %3.1f +/- %2.1f \t\t [%3.1f-%3.1f] \t ms \n',easystats(ICI)) ;
fprintf(' RANGE \t\t\t\t %3.1f +/- %2.1f \t\t [%3.1f-%3.1f] \t m \n',easystats(RANGE)) ;
fprintf(' N \t\t\t\t\t\t %3.0f \n',length(ASLpp)) ;

fprintf(['\n Run pistonfit.m with input ' datafile ' to get directional characteristics of sonar beam \n\n'])

% Find highest back-calculated signal and return as model on-axis click
kmax    = find(ASLpp==max(ASLpp)) ;
sig       = SIGNAL(kmax,:) ;
fs      = fs(kmax) ;
