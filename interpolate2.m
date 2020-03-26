function [o]=interpolate2(in,N)
L=length(in);
IN=fft(in);
IN(round(length(IN)/2)+1:floor(N*L))=0;
if N<1,IN(round(N*L/2)+1:L)=0;end;
IN=N*IN(1:floor(N*L));
IN(1)=IN(1)/2;
o=2*real(ifft(IN));
