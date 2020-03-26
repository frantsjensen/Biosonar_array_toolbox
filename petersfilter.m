% ref: 1MHz
zf = 0.150 ;  % zero frequency = 0.15 MHz
pf = 0.4 ;   % pole frequency = 0.4 MHz
Fs=1;

A = pf/zf*[1 -2*pi*zf] ;      % numerator
B = [1 2*pi*pf] ;            % denominator

[h,w]=freqs(A,B,2*pi*linspace(0.01,1,1024));

semilogx(w/2/pi*1000,20*log10(abs(h))),grid

% help bilinear

% 2nd order prototype
% A = [1 -2*2*pi*zf/q (2*pi*zf)^2]

[NUMd,DENd] = BILINEAR(A,B,Fs);
load Frants
H=freqz(NUMd,DENd,frequency,1e6);
frequency=[frequency]*1000;
voltage=20*log10([voltage]/965);
plot(frequency,20*log10(abs(H))+voltage)