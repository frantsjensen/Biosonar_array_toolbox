function [CMAT] = channelcolors(N)
% Quick script to define colors for some click parameter tools
%
% F. H. Jensen, 2013 (frants.jensen@gmail.com)

% Define colors for plot - only 8 colors defined per default
CMAT = [238 045 046 ;... % Red
        024 089 169 ;... % Blue
        000 139 071 ;... % Green
        243 125 035 ;... % Orange
        105 042 147 ;... % Purple
        162 029 032 ;... % Rust
        179 056 147 ;... % Pink/purple
        001 001 001]/255;% Black

% Let's make color matrix really big in case someone decides to use click
% parameter toolbox for arrays with more than 8 receivers
CMAT = [CMAT ; CMAT ; CMAT; CMAT; CMAT; CMAT ];
    
% Return only the amount of different colors needed
CMAT = CMAT(1:N,:);