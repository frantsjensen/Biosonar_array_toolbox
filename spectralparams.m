function [Fparams Fpsd] = spectralparams(sig,fs,nfft)

% Calculate spectral parameters for a click-like (short-duration) signal
%
% Input: 
%   sig:    signal waveform to be evaluated
%   fs :    sample rate of waveform
%   nfft:   Desired FFT size of discrete fourier transform
%
% Signal should be a window of L samples centred around the peak of the signal.
% Window length L should be a whole number to the power of 2 ([2, 4, 8, 16, 32, 64, etc].
% Perform desired frequency filtering and windowing on signal before
% running spectralparams - only a default rectangular window is included.
%
% If fft size is larger than length of signal, the fft will be interpolated
% with a high-quality Sinc interpolation (non-linear) before parameters are
% estimated. The window length L will still define the spectral resolution
% so that spectral resolution df is given by:
%   df = (fs/2) / (L/2) = fs/L 
% Interpolated spectral resolution (interval between frequency points) is
% given by df' = fs / nfft
%
% Output:
%   Fparams: Vector of the following parameters (all in kHz):
%       Peak frequency (Fp), centroid frequency (fc), Centralized RMS bandwidth
%       (Bwrms), -3dB (BW-3dB) and -10 dB (BW-10dB) bandwidth.
%       Fparams = [Fp Fc BWrms BW-3dB BW-10dB]
%   Fpsd : Matrix of frequency (first column, kHz) and power spectral
%       density (second column). Second column needs to be corrected for
%       peak clip level
%
% F. H. Jensen, 2013 (frants.jensen@gmail.com)
% Based on Madsen and Wahlberg 2007 equations

% Check that signal is vector
if length(sig(:))~=length(sig)
    error(' Input signal needs to be vector')
end

% First, define fft size of discrete fourier transform
% If N is larger than length of signal, spectrum will be 
% interpolated (Sinc interpolation)
L = length(sig) ;
N = nfft ; 
interpolation = N/L ;

% Apply rectangular window to make sure end points are 0
sig = sig(:).*rectwin(L) ;

% Calculate amplitude spectrum and power spectral density using discrete fft
fsig     = abs(fft(sig,N));                                 % Calculate amplitude spectrum
fsig     = fsig(1:length(fsig)/2);                          % Drop second half of amplitude spectrum
fpsd     = 2 * (1/fs/(N).*abs(fsig).^2) ;                   % Power spectral density estimate
fsig     = fsig/max(fsig);                                  % Normalize amplitude spectrum
f        = (1/1000)*(fs/2)*[0:length(fsig)-1]'/(length(fsig)); % Frequency vector, kHz

% Calculate frequency parameters
fp    = f(find(fsig==max(fsig)));                        % Peak frequency in kHz
fc    = sum(f.*fsig.^2)/sum(fsig.^2);                    % Centroid frequency, kHz
bwrms    = sqrt(sum((f-fc).^2.*fsig.^2 ) / sum(fsig.^2)); % RMS centralized bandwidth, kHz

f_ind    = find(fsig.^2>10^(-3/10));                        % Find index of power spectrum within -3dB of peak power
bw3db    = (f(max(f_ind))-f(min(f_ind)));                   % -3 dB BW in kHz
f_ind    = find(fsig.^2>10^(-10/10));                       % Find index of power spectrum within -10dB of peak power
bw10db   = (f(max(f_ind))-f(min(f_ind)));

% Accumulate frequency parameters in one vector
Fparams = [fp fc bwrms bw3db bw10db] ;

% Build matrix of frequency (kHz, column 1) and 
% power spectral density (column 2)
Fpsd = [f(:) fpsd(:)] ;
