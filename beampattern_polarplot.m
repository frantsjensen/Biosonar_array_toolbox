function [] = CreateStenellaPP()

cd 'c:\matlab_scripts\AGC\2011 Stenella';
D=load('bootstrap_stenella_logfit.mat');

% Polar plot settings
dB_range = 30;
span = [-45 45];
orientation = 90;
ppstyle = ':';

figure(1), clf, hold on, box on
%polar_plot_line_LJ(30,[-45 45],90,D.pist_model_lower,D.piston_angles','y',[0.6 0.6 0.6],2)
%polar_plot_line_LJ(30,[-45 45],90,D.pist_model_upper,D.piston_angles','y',[0.6 0.6 0.6],2)
%polar_plot_line_LJ(30,[-45 45],90,D.pist_model_db,D.piston_angles','y',[0.2 0.2 0.2],2.5)
%set(gca,'xlim',[0 30],'XTick',[],'YLim',0.5*[-30 30],'YTick',[-30:6:30])


% Create polar frame
makepolarframe (dB_range,span,orientation,ppstyle) ;

% Now add data to polar plot
h1 = addpolarline (D.pist_model_lower,D.piston_angles',dB_range,'y',90,[0.5 0.5 0.5],1.5,'--');
h2 = addpolarline (D.pist_model_upper,D.piston_angles',dB_range,'y',90,[0.5 0.5 0.5],1.5,'--');
h3 = addpolarline (D.pist_model_db,D.piston_angles',dB_range,'y',90,[0 0 0],1.5);

% Adjust axis settings (axes need to be equal)
set(gca,'xlim',[-8 40],'XTick',[],'YLim',[-24 24],'YTick',[-30:6:30])
set(gca,'YTick',[]), set(gca, 'visible', 'off') ; 
axis equal

% Set paperunits to inches
set(gcf,'PaperUnits','inches')

% Find current papersize
papersize = get(gcf, 'PaperSize');

% Adjust desired width and height
width = 1.5*3.375;    % JASA 1-column width: 3.3750 inch
height = 1.5*3.3750;          

left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);

%print -dtiff -r300 Jensen_stenella_Fig4_alt




function [hpp] = makepolarframe (dB_range,span,orientation,ppstyle)
%dB_range = magnitude range eg. 30 creates a plot spanning 0 to -30 dB.
%span = angular span of plot in degrees eg. [-90 90] creates standard 180 degree upwards plot
%orientation = angular rotation of the plot in degrees
%ppstyle = line style of polar plot framework

if nargin<4,
    ppstyle = 'k';
end

if nargin<3,
    orientation = 0;
end

% Create polar plot framework
Angle_range = (span(1):15:span(2))+orientation;
[y_data,x_data] = pol2cart(Angle_range*pi/180,dB_range + 3);
alpha = -dB_range:6:0;
theta = Angle_range(1)+180:0.1:Angle_range(end)+180;
for m = 1:length(alpha)
        [opp(:,m) adj(:,m)] = pol2cart(theta*pi/180,alpha(m));
end
% add polar plot framework to figure;
for n = 1:length(adj(1,:))
plot(adj(:,n),opp(:,n),'color',[0.5 0.5 0.5],'LineStyle',ppstyle)
hold on
end
for n = 1:length(x_data)
    hpp(n)=plot([0 x_data(n)],[0 y_data(n)],'color',[0.5 0.5 0.5],'LineStyle',ppstyle);
end
plot([0 x_data(1)],[0 y_data(1)],'k')
plot([0 x_data(end)],[0 y_data(end)],'k')





function [h] = addpolarline (magnitude,angle,dB_range,mirror,orientation,col,width,style)
%Create lines for data within polar plot
%magnitude = relative amplitude measurements 
%angle = angle measurements
%mirror = if 'y', the results are mirrored around 0 degrees, eles write 'n'
%orientation = orientation of polar plot
%dB_range = dB range of polar plot
%col = color of line
%width = linewidth
%style = line style

if nargin<8,
    style = '-';
end

if nargin<7,
    width = 2;
end

if nargin<6,
    col = 'k';
end

if nargin<5
    orientation = 0;
end

if nargin<4,
    mirror='y';
end

angle=angle(:); magnitude=magnitude(:);

if mirror == 'y',
    angle = [-flipud(angle) ; NaN ; angle];
    magnitude = [flipud(magnitude) ; NaN ; magnitude];
end

[y_points x_points] = pol2cart((angle+orientation)*pi/180,magnitude+dB_range);
h=plot(x_points,y_points,'color',col,'linewidth',width,'LineStyle',style);
