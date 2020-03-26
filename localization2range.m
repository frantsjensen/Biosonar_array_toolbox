function range = localization2range(LOC,R)
% Modular script to estimate range between source localization in 2D or 3D
% (functionality needs to be added depending on receiver geometry) and a
% linear or 2D array
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
    % Calculate range (meters) for 1D vertical array from 2D localization 
    % (horizontal range and depth)
    range = ( sqrt( LOC(1)^2 + (LOC(2)-R).^2) )';
    
elseif length(LOC)==3, 
    %To make this, R (receiver coordinates) probably needs to be 2d, 
    % meaning a change to toolbox_settings.m and an update to
    % localize_click.m
    error('Error in localization2range - needs to be reconfigured for 3D localization');
end