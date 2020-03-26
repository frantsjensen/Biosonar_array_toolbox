% Calculates and returns the p-p and rms received source levels (in dB re 1uPa), 
% the time window (in u-seconds) and the Energy flux density (in uPa^2*s)
%
% Input requires the signal, the sampling frequency and the absolute value
% clip level (the received sound pressure level at the peak value of 1,
% given in dB re 1uPa) as well as a request for the desired method (either 1 or 2)
% and returns the peak-to-peak and RMS received level corrected 
% for amplification as well as the window in u-seconds and the energy flux
% of the pulse.
%
% Function takes the following form:
% [ peak_peak , RMS  ,time_window , Energy_flux ] = RL (sig , fs , abs_clip_level , method)
%
% Method indicates whether RMS, Tau and Energy flux are based on
% a) a -10 dB amplitude criterion (for method = 1)
% b) a 95% Energy criterion (for method = 2)
%
% Last edited 2013 by F. H. Jensen (frants.jensen@gmail.com)
% Based on Madsen and Wahlberg 2007 equations

function [ peak_peak , RMS , time_window , Energy_flux ] = RL (sig , fs , clip_level , method )

if nargin<4,
    method = 1; % Default method is -10 dB amplitude criterion
end

if nargin<3
    clip_level = 0;
end

%Testing for length of signal
if length(sig) <= 10 * 1e-6 * fs, disp('Warning - signal shorter than 10 us'),end

%Received source levels in dB re 1µPa
peak_peak = clip_level + 20*log10( max(sig) - min(sig) );

if method == 1
    %Choosing time window - -10dB amplitude criterion
	sig10       = interp ( sig , 10);                       %Interpolates 10x
	hill10      = sqrt(sig10.^2 + imag(hilbert(sig10)).^2);
	I10         = find ( hill10 >= max(hill10)/(10^(10/20)));   %Find -10 dB energy indexes [correct]
	time_window = 1e6 * length (I10) / (10*fs);             %Length of (-10 dB) signal (u-seconds)

    %Calculating actual RMS sound pressure and Energy flux   
    RMS = clip_level + 20*log10( sqrt ( mean (sig10 (I10).^2)) );
    Energy_flux = RMS  + 10*log10 ( 1e-6 * time_window );

elseif method == 2
    %Choosing time window - 95% Energy criterion
    criterion    = 0.95;                    %Sets the energy criterion
    
    % Interpolate signal
	sig10        = interp ( sig , 10);
    
    % Find 95% percentile time window
    total_energy = sig10'*sig10;            %Calculates the summed energy
	energy       = sig10.^2;                %Instantaneous energy
	for i=1:length(sig10),
        if sum ( energy(1:i) ) < total_energy * (1-criterion)/2, e_low=i; end
        if sum ( energy(1:i) ) < total_energy * (1 - (1-criterion)/2), e_high=i; end
	end    
	time_window  = 1e6 * ( e_high-e_low ) / (10*fs) ;
	sigfinal     = sig10 (e_low : e_high) ;

%     figure, hold on, 
%     plot([1:length(energy)]*1e6/(10*fs),cumsum(energy),'k'), 
%     plot([1 length(energy)]*1e6/(10*fs),[1 1]*total_energy,'b')
%     plot([1 e_low e_low]*1e6/(10*fs), [1 1 0]*total_energy*(1-criterion)/2,'r')
%     plot([e_high e_high length(energy)]*1e6/(10*fs), [0 1 1]*total_energy*(1-(1-criterion)/2),'r')
%     plot([1:length(energy)]*1e6/(10*fs),0.5*total_energy.*(sigtemp+1),'y')
    
    %Calculating corrected RMS sound pressure and Energy flux   
	RMS = clip_level + 20*log10( sqrt ( mean ( sigfinal.^2)) );
	Energy_flux = RMS  + 10*log10 ( 1e-6 * time_window );
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Script