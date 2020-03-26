function [pt,varargout]=lagrange3pt(varargin)
%pt=lagrpt([x0 x1 x2],[y0 y1 y2])
%returns the peak amplitude (if a peak is present between x0 and x2) 
%of the lagrange polynomium interpolating the three number pairs
% Kristian Beedholm, Aarhus University, 2007

if nargin==6
    x0=varargin{1};x1=varargin{2};x2=varargin{3};
    y0=varargin{4};y1=varargin{5};y2=varargin{6};
elseif nargin==2;
    xar=varargin{1};
    yar=varargin{2};
    x0=xar(1);x1=xar(2);x2=xar(3);
    y0=yar(1);y1=yar(2);y2=yar(3);
else
    error('bad input format!')
end;
    
taeller=(y0*x2^2-y0*x1^2-y1*x2^2+y1*x0^2+y2*x1^2-y2*x0^2);
naevner=y0*x2-y0*x1-y1*x2+y1*x0-y2*x0+y2*x1;
pt=0.5*taeller/naevner;

pk=((pt-x1)*(pt-x2)/((x0-x1)*(x0-x2)))*y0+...
   ((pt-x0)*(pt-x2)/((x1-x0)*(x1-x2)))*y1+...
   ((pt-x0)*(pt-x1)/((x2-x0)*(x2-x1)))*y2;


if nargout>1
varargout={pk};
end;