function pist_model_db = simulatebeampattern(piston_size,piston_angles,s,nsr)
% pist_model_db = simulatebeampattern(piston_size,piston_angles,s,nsr)
% Calculate beam pattern based on circular piston model
% 
% Input:
% - piston_size     Piston diameter of circular piston (in cm)
% - piston_angles   Angles to calculate beampattern for
% - s               interpolated on-axis signal to convolve
% - nsr             sampling rate of interpolated signal
%
% output will be power in dB relative to on-axis power
% F. H. Jensen, 2013 (frants.jensen@gmail.com)

for j=1:length(piston_angles) 
    pist_model(j)=energy(conv(imppist(piston_angles(j),piston_size/100,1.5,nsr),s));
end
%pist_model=pist_model/max(pist_model); %normalize to max
%pist_model_db=smooth(piston_angles,10*log10(pist_model),0.1,'rloess');

pist_model=smooth(piston_angles,pist_model,0.1,'rloess');
pist_model=pist_model/max(pist_model); %normalize to max
pist_model_db=10*log10(pist_model);
