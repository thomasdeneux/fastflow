function [Y, CS] = fast_loaddata_ref(Yname,fname,varargin)
%---
% function Y = fast_loaddata_ref('Y',fname,'headers')
% function [Y CS] = fast_loaddata_ref('Y',fname,[frames[,sux,suby]|[,indices]][,'normalize'])
%---
% ne doit pas recr???er de l'espace m???moire - trop chiant ce faux passage par
% r???f???rence
%
% See also fast_loaddata

if nargin==0
    disp('usage:')
    help fast_loaddata_ref
    return
end

Y = evalin('caller',Yname);
evalin('caller',['clear ' Yname])
narg = nargin-1;

% change working directory
swd = pwd;
%switch fname(1:2)
 %   case 'TC'
        catflag = false;
        datatype = 0;
 %       switch fname(20)
 %           case '1'
 %               cd /storm5/thomas/Data/0503_Marseille/TOM03oct05
 %           case '2'
 %               cd /storm5/thomas/Data/0503_Marseille/TOM03oct05
 %           case '3'
 %               cd /storm5/thomas/Data/0503_Marseille/TOM03oct05
 %       end
 %   case 'fl'
 %       catflag = true;
 %       datatype = 2;
 %       cd C:\Users\Thomas\WORK\IRMf\0503_Marseille\data\catflow\
% end

% get headers
[nfrms,xs,ys,hsize,nstims]=get_head(fname,'nframesperstim,framewidth,frameheight,lenheader,nstimuli');

% headers only
if narg==2 && any(strfind(varargin{1},'header'))
    Y = struct('nfrms',nfrms,'xs',xs,'ys',ys,'hsize',hsize,'nstims',nstims);
    return
end

% Input
% frames
if narg>1 && ~ischar(varargin{1})
    frames = varargin{1};
    nt = length(frames);
else
    frames = 1:nfrms;
end
nt = length(frames);
% sub-window
if narg>2 && ~ischar(varargin{2})
    if narg>3 && ~ischar(varargin{3})
        subx = varargin{2};
        suby = varargin{3};
        [yy xx] = meshgrid(suby,subx);
        indices = xx + xs*(yy-1);
    else
        indices = varargin{2};
    end
else
    subx = 1:xs; suby = 1:ys;
    indices = 1:xs*ys;
end

% load data
Y = loadsum2_ref('Y',fname,frames-1,0,datatype,xs,ys,hsize,indices);

% normalize
if narg>1 && any(strfind(varargin{end},'norm'))
    m = mean(Y,2);
    for i=1:size(Y,1)
        Y(i,:) = Y(i,:)./m(i);
    end
    if nargout>1 && exist('subx','var')
        CS = reshape(m,length(subx),length(suby));
    end
end
if exist('subx','var')
    Y = reshape(Y,length(subx),length(suby),nt);
end

% back to old working directory
cd(swd)
