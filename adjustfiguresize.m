function [] = adjustfiguresize(width,height,units)
% Quickly adjust figure size for high-quality print
%
% adjustfiguresize(width,height)
% sets the current figure window to a given size (width, height) in inches
% most scientific journals (such as jasa) has a 1-column width of 3.375 inches
% for good aspect ratio, use either square figure (height = width)
% or golden ratio (height = 1.618 * width or height = width / 1.618)
%
% MatLab does not render well to tiff or other rasterized formats
% To circumvent this, print to double required size and downsample in image
% processing program (Adobe Photoshop or similar)
%
% This script does not adjust font type or size, tick marks, or other plot options
% MatLab default is to change tick marks when expanding figure, which does
% not translate well into a rasterized image. These can be set manually by
% changing axis properties ('XTick' and 'YTick'). Manual tick labelling can
% also be done by changing axis properties 'XTickLabel' and 'YTickLabel'.
%
% Example code for changing other figure formats here:
% To change current axis fonts and axis labels fonts:
% - set(gca,'FontName','Helvetica','FontSize',12)
% - handle_labels = get(gca,{'XLabel','YLabel'});
% - set([handle_labels{:,:}],'FontName','Helvetica','FontSize',14)
% To change ticks:
% - set(gca,'XTick',xticks,'YTick',yticks)
% To make box around plot:
% - box on
%
% F. H. Jensen, 2013 (frants.jensen@gmail.com)

if nargin<1
    width = 6.750 ; % Default width 2-column format
end

if nargin<2
    height = width/1.618 ; % Default aspect is landscape, golden ratio
end

if nargin<3
    units = 'inches' ; % Default unit is inches
end

% Set paperunits to inches
set(gcf,'PaperUnits',units)

% Find current papersize
papersize = get(gcf, 'PaperSize'); 

left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);