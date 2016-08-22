function params = fast_params(varargin)
% function params = fast_params(['tensor'][,interpstep,hight,highx,lowt,lowx,sigmat,sigmax])
% function params = fast_params('fft'[,interpstep,hight,highx,lowt,lowx,windowwidth,windowcover,sigmaangle,searchangle])
%---
% returns structure with parameters (to use in fast_dotensor or fast_dofft)

if nargin==0, help fast_params , return, end 

if nargin==0
    mode='tensor';
    idx=0;
else
    mode=varargin{1};
    idx=1;
end

switch mode
    case 'tensor'
        fields={'interpstep','hight','highx','lowt','lowx','sigmat','sigmax'};
        defaults=[1 0 30 2.5 2.5 30 5];
    case 'fft'
        fields={'interpstep','hight','highx','lowt','lowx','windowwidth','windowcover','sigmaangle','searchangle'};
        defaults=[1 50 50 0 0 25 1 25 1/8];
end
nfields = length(fields);

if nargin==idx
    for i=1:nfields
        params.(fields{i})=defaults(i);
    end
else
    if nargin-idx~=nfields, error('wrong number of parameters (%i instead of %i)',nargin-idx,nfields), end
    for i=1:nfields
        params.(fields{i})=varargin{i+idx};
    end
end
