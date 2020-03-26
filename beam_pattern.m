function pist_model_db = beam_pattern(piston_size,piston_angles,s,nsr)

% Calculate beam pattern based on circular piston model
% Piston diameter (in cm), s is the interpolated on-axis signal to convolve, 
% nsr the interpolated sampling rate in Hz
% output will be relative power in dB

try 
    SETTINGS = toolbox_settings();
    soundspeed = SETTINGS.soundspeed;
catch
    soundspeed = 1524 ; % Default value
end

for j=1:length(piston_angles) 
    pist_model(j)=energy(conv(imppist(piston_angles(j),piston_size/100,soundspeed,nsr),s));
end

pist_model=smooth(piston_angles,pist_model,0.1,'rloess');
pist_model=pist_model/max(pist_model); %normalize to max
pist_model_db=10*log10(pist_model);
