% Change code here to locate correct file with compiled click parameters
% Go through plot functions gradually (mark lines and run with F9)
% some axis limits MUST be adjusted to better fit your data
cd 'c:/Matlab_scripts/AGC/2011 Stenella/'
load stenella_revised.mat


%%%%%%%%%%%%%%% FIGURE A %%%%%%%%%%%%%%%

% Description:
% Plot the interpolated waveform of 4-5 high-quality on-axis signals

% Isolate 4 clicks with highest backcalculated ASLpp
[sortasl k]=sort(ASLpp);
k=find(ASLpp>=sortasl(end-5));

% Plot signal interpolated x100 and normalized
figure(1), clf, 
subplot('position',[0.15 0.63 0.8 0.32]), hold on
for i=1:5,
    sig = interp(SIGNAL(k(i),:),100);
    hilb = abs(hilbert(sig));
    peak = find(hilb==max(hilb));
    h1(i)=plot( ([1:length(sig)]-peak)/(fs*100/1e6) , sig/max(abs(sig)) , 'k');
    % MAYBE CALCULATE ABS AMPLITUDE INSTEAD
    set(h1(i),'Color',(1-0.25*i)*[1 1 1])
end
    
% Adjust plot
box on,
set(h1,'LineWidth',2)
set(gca,'Xlim',[-30 30])
xlabel('Time, µs','FontSize',14,'FontName','helvetica')
ylabel('Normalized amplitude','FontSize',14,'FontName','helvetica')



%%%%%%%%%%%%%% FIGURE B %%%%%%%%%%%%%%%%

% Plot normalized spectra and average spectrum
figure(1), 
subplot('position',[0.15 0.1 0.8 0.45]), hold on
freq=FREQ*1000;
for i=1:size(SPEC,1), sp=SPEC(i,:)/max(SPEC(i,:); h1=plot(freq/1000,10*log10(sp),'k','Color',0.4*[1 1 1]); end
h2=plot(freq/1000,10*log10(mean(SPEC)),'k','Color',[0.0 0.0 0.1],'LineWidth',4);
box on, set(gca,'YLim',[-25 1],'XLim',[0 200])
%legend([h2],'Mean power spectrum')
xlabel('Frequency, kHz','FontSize',14)
ylabel('Relative intensity, dB','FontSize',14)

% Add a box plot of Centroid frequency estimates to the figure (this is
% unnecessary, but see how it can be used in Wahlberg et al to show
% frequency difference between two species)
[whiskers,whiskerends,mybox,medianline]=boxutil_fhj(Fc,-20,4);
hpoints=plot(Fc,-21.5*ones(length(Fc),1)+3*rand(length(Fc),1),'ko','MarkerEdgeColor',[0.7 0.7 0.8],'MarkerFaceColor',[0.7 0.7 0.8],'MarkerSize',2);
plot(whiskerends(2,:),whiskerends(1,:),'k-','LineWidth',2);
plot(whiskers(2,:),whiskers(1,:),'k--','LineWidth',2);
plot(mybox(2,:),mybox(1,:),'k','LineWidth',2);
plot(medianline(2,:),medianline(1,:),'r','LineWidth',2,'Color',[0.5 0 0]);
%text(85,0.5,'A: Power spectra','HorizontalAlignment','Center','VerticalAlignment','Bottom')
%text(85,-22,'B: Centroid frequency','HorizontalAlignment','Center','VerticalAlignment','Bottom')
set(gca,'fontsize',12)


% Set paperunits to inches
set(gcf,'PaperUnits','inches')

% Find current papersize
papersize = get(gcf, 'PaperSize');

% Adjust desired width and height
width = 2*3.375;    % JASA 1-column width: 3.3750 inch
height = 2*3.3750;          

left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);

adjustfigurefont('Helvetica',12,14)

%print -dtiff -r300 myfigurename


%%%%%%%%%%%%%% FIGURE C %%%%%%%%%%%%%%%%

% RANGE vs SL
figure(2), clf,
subplot('position',[0.15 0.58 0.8 0.36]), hold on, box on
[P,S]=polyfit(log10(RANGE),ASLrms_amp,1);
Xpred=log10([0.1:0.1:1 2:50]);
[Y,delta]=polyval(P,Xpred,S);
fill([10.^(Xpred) fliplr(10.^Xpred)],[Y-delta fliplr(Y+delta)],[0.9 0.9 0.9],'LineStyle','none','FaceColor',[0.9 0.9 0.9],'FaceAlpha',1)
plot(10.^Xpred,Y,'k','LineWidth',2,'Color',[0.5 0.5 0.5])
plot(RANGE,ASLrms_amp,'ks','MarkerFaceColor',[0 0 0],'MarkerSize',5)
set(gca,'FontSize',12,'layer','top')
xlabel('Range, m','FontSize',14)
ylabel('ASL, dB re. 1 µPa rms','FontSize',14)
text(3,207,'A','FontSize',16,'HorizontalAlignment','Right')
disp(['ASL = ' num2str(P(1)) ' * log10(R) + ' num2str((P(2)))])

% Fc versus SL
subplot('position',[0.15 0.10 0.8 0.36]), hold on, box on
[P,S]=polyfit(ASLrms_amp,Fc,1);
Xpred=[185:210];
[Y,delta]=polyval(P,Xpred,S);
fill([Xpred fliplr(Xpred)],[Y-delta fliplr(Y+delta)],[0.9 0.9 0.9],'LineStyle','none','FaceColor',[0.9 0.9 0.9],'FaceAlpha',1)
plot([185 210],P(1)*[185 210]+P(2),'k','LineWidth',2,'Color',[0.5 0.5 0.5])
plot(ASLrms_amp,Fc,'ks','MarkerFaceColor',[0 0 0],'MarkerSize',5)
disp(['Fc = ' num2str(P(1)) ' * ASL ' num2str(P(2))])

set(gca,'FontSize',12,'layer','top','YLim',[60 100])
xlabel('ASL, dB re. 1 µPa rms','FontSize',14)
ylabel('Fc, kHz','FontSize',14)
text(186.4,96,'B','FontSize',16,'HorizontalAlignment','Right')


% Set paperunits to inches
set(gcf,'PaperUnits','inches')

% Find current papersize
papersize = get(gcf, 'PaperSize');

% Adjust desired width and height
width = 1.5*3.375;    % JASA 1-column width: 3.3750 inch
height = 2*3.3750;          

left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);

%print -dtiff -r300 myfigurename2


%%%%%%%%%%%%%% FIGURE D %%%%%%%%%%%%%%%%

% Take only closest signals
k=find(RANGE<20);

%D=load('stenella_pistonfit.mat');
D=load('bootstrap_stenella_logfit.mat');
figure(3), clf, hold on, box on,
h2 = plot(abs(ANGLES(k,:)),RELASL(k,:),'ks','MarkerFaceColor','k','MarkerSize',3);
h2(2) = plot(D.piston_angles,D.pist_model_best,'k','LineWidth',3,'color',[0.3 0.3 0.3]);
h2(3) = plot(D.piston_angles,D.pist_model_lower,'k--','LineWidth',2,'color',[0.3 0.3 0.3]);
h2(4) = plot(D.piston_angles,D.pist_model_upper,'k--','LineWidth',2,'color',[0.3 0.3 0.3]);

%polar(abs(2*pi*ANGLES(k,:)/360),40+ASLall(k,:)-[AMP(k);AMP(k);AMP(k);AMP(k);AMP(k);AMP(k)]','kx')
set(gca,'fontsize',12,'ylim',[-25 1],'xlim',[0 25])
xlabel('Angle of incidence (degrees vertical)','FontSize',14)
ylabel('ASL (dB re. 1 µPa pp)','FontSize',14)
legend(h2(1:3),'Array data','Piston fit','Piston fit bootstrap CI','Location','SouthWest')
legend boxoff
disp(['Clicks within 20m: ' num2str(length(k))])


% Set paperunits to inches
set(gcf,'PaperUnits','inches')

% Find current papersize
papersize = get(gcf, 'PaperSize');

% Adjust desired width and height
width = 1.5*3.375;    % JASA 1-column width: 3.3750 inch
height = 1.2*3.3750;          

left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);

% print -dtiff -r300 myfigurename3






