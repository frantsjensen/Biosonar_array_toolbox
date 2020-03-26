function [] = adjustfigurefont(fontname,fontsize,fontsizelabels)
% Quickly adjust figure font for high-quality print
% adjustfigurefont(fontname,fontsize,fontsizelabels)
% will change font used on all axes and all axis labels in current figure
% Default values are
% fontname = 'Helvetica' ;
% fontsize = 12 ;
% fontsizelabels = fontsize + 2 ;

if nargin<1,
    fontname = 'Helvetica' ; % Default font Helvetica
end

if nargin<2,
    fontsize = 12 ; % Default font size 12
end

if nargin<3
    fontsizelabels = fontsize + 2 ; % Default font size for labels 14
end

% Find handles for all axes within current figure
handle_axes = get(gcf,'children') ;

% Change font on those axes
set(handle_axes,'FontName',fontname,'FontSize',fontsize)

% Find handles for all axis labels associated with those axes
handle_labels = get(handle_axes,{'XLabel','YLabel'});

% Change font on those labels
set([handle_labels{:,:}],'FontName',fontname,'FontSize',fontsizelabels)
