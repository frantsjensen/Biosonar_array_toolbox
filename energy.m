function [o]=energy(varargin);
if nargin==0
    s=coolcopy;
else
    s=varargin{1};
end;

o=sum(s.^2);