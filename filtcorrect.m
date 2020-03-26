%Filter to correct for freq response in  NUK 4 ch filter&amp box using the
%200 kHz LP filter. See peterfilter.m for how coef are derived.

function [sigf] = filtcorrect (sig,fs);

%file='c:/temp/coolclipboard.wav';  %load signal
%[sig,fs,nbit]=wavread(file);
x=interp(sig,2);  %interpolate to 1MHz to allow for second LP filter to operate

NUMd = [0.6248 -1.7386];

DENd = [1.0000 0.1137];  %filter coefficients
    
y=filter(NUMd,DENd,x); %apply filter

sigf=decdc(y,2);  %downsample to 500k

%wavwrite(-sigf,fs,nbit,'c:/signalf.wav');  %write wave file