% Fill in data from array_calibration.m for each calibration range (file)
% For each calibration range, run through array_calibration.m
% Find rms error for each range, then convert to relative rms error
%
% F. H. Jensen 2013 (frants.jensen@gmail.com)

% Range and error estimates
dist = [5 10 20 30 40]; % True range
range = [5.02 10.0492 21.7538 34.36 35.25]; % Mean range estimate
st_dev = [0.08 0.197 0.945 2.48 5.35];      % Range estimate standard deviation
RMS_error_percent = [0.016 0.02 0.099 0.167 0.178]*100; % Range estimate rms error
TL_RMS_error = -20*log10(1-RMS_error_percent/100);     % Transmission loss rms error

subplot('position',[0.15 0.15 0.8 0.35]),
h3 = plot(dist,TL_RMS_error,'k.','MarkerSize',16);
set(gca,'FontSize',12, 'FontName' , 'Helvetica','xlim',[0 45],'Ylim',[0 3],'Ytick',[0:1:3])
lab3 = xlabel(char('True range (m)'), 'FontName' , 'Helvetica' , 'FontSize' , 12);
temp = get(lab3,'position');temp(2)=-1.20;
%set (lab3,'position',temp)
grid on

lab3 = ylabel(char('Transmission loss','  RMS Error (dB)'), 'FontName' , 'Helvetica' , 'FontSize' , 12);
temp = get(lab3,'position');temp(1)=-2.8;
set (lab3,'position',temp)

subplot('position',[0.15 0.55 0.8 0.35]), hold on,
plot([0 40],[0 40],'k--','LineWidth',1)
h1 = errorbar(dist,range,st_dev,st_dev,'k.');
hold off

set(h1,'MarkerSize',16)
set(gca,'FontSize',12, 'FontName' , 'Helvetica','xlim',[0 45],'Ylim',[0 50],'Xticklabel',[], 'Ytick',[0:10:40], 'box', 'on')
lab1 = ylabel(char('Localization', 'mean range (m)'), 'FontName' , 'Helvetica' , 'FontSize' , 12);
temp = get(lab1,'position');temp(1)=-4;
grid on

text(dist,(max(dist)+5)*ones(1,length(dist)),['N=600';'N=695';'N=695';'N=634';'N=196'],'FontName','Helvetica','FontWeight','Bold','HorizontalAlignment','center','verticalAlignment','middle','backgroundcolor','w')

adjustfiguresize(2*3.375,1.618*2*3.375)

% Now print to 300 dpi tiff then change format to 1-column width (downsample):
print -dtiff -r300 myarraycalibration