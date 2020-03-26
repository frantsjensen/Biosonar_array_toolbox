% Apply shading to figure
%
% []=shade(x,y,color)
function []=shade(x,y,color);
fill(x,y,color,'LineStyle','none','FaceColor',color,'FaceAlpha',0.5)