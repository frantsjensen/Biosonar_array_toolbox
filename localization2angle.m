function theta = localization2angle(LOC,R,ASLrms)
% Modular script to estimate angle of incidence between source localization in 2D or 3D
% and a linear or 2D array ( for 2d array needs to be added depending on receiver geometry)
%
% Input:
% LOC       Column vector with source localization
% R         Column vector or matrix with receiver coordinates
%
% 1D array: LOC will be a 2-element column
% vector consisting of vertical range and depth
% R will be an Nx1 column vector consisting of depth of each hydrophone
%
% 2D array (needs to be rewritten): LOC will be a 3-element column vector
% (likely [z,y,x]' coordinates where z is the normal vector to the plane of
% the hydrophones, y is depth, and x is the last dimension of the array 
% - but this will be defined once localization routine is adjusted) 
%
% F. H. Jensen, 2014 (frants.jensen@gmail.com)

if length(LOC)==2,
    % Calculate angle of arrival for 1d array
    theta = zeros(1,length(R));
    for j=1:length(R), 
        theta(j) = 0.5*pi + sign(LOC(2)-R(j)) * atan (abs(R(j)-LOC(2)) / LOC(1) );  %Radians
    end  
    theta=(180/pi)*theta; % Angle in degrees
    
    
elseif length(LOC)==3, 
    %To make this, R (receiver coordinates) probably needs to be 2d, 
    % meaning a change to toolbox_settings.m and an update to
    % localize_click.m
    error('Error in localization2angle - needs to be reconfigured for 3D localization');
end