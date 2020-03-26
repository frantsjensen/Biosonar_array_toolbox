function [eh]=imppist(dangle,D,c,sr)
% Returns impulse response of circular piston
% Input:
% -  dangle  off-axis angle in degrees
% -  D       piston diameter (m)
% -  c       sound speed of medium (m/s)
% -  sr      sample rate (Hz)
% Output:
% -  eh      time-domain impulse response of piston
%
% Instructions:
% Convolve on-axis signal (possibly upsampled) with impulse response to get
% expected off-axis signal
%
% Script by K. Beedholm, Aarhus University

alpha=dangle/180*pi; %rad

if dangle>0 
   L=sr/c*(sin(alpha)*D);
   if floor(L)>1
      fakts=sin(acos(2*([0:round(L)-1])/L-1));
      eh=fakts(1:length(fakts))/sum(fakts); %normalize
   else
      eh=1;
   end;
else
   eh=1;
end;
