function [Y, CS] = fast_loaddata_single(fname,varargin)
% same as fast_loaddata, but returns arrays of singles
%---
% function Y = fast_loaddata(fname,'headers')
% function [Y CS] = fast_loaddata(fname,[frames[,sux,suby]|[,indices]][,'normalize'])
%
% See also fast_loaddata_ref

if nargin==0
    disp('usage:')
    help fast_loaddata
    Y = [];
    return
end

disp('change working directory disabled - assumed datype = 0')
datatype = 0;
swd = pwd;
catflag = false;
% % change working directory
% swd = pwd;
% datatype = 0;
% switch fname(1:2)
%     case 'TC'
%         catflag = false;
%         datatype = 0;
%         switch fname(20)
%             case '1'
%                 cd /storm5/thomas/Data/0503_Marseille/TOM03oct05
%             case '2'
%                 cd /storm5/thomas/Data/0503_Marseille/TOM03oct05
%             case '3'
%                 cd /storm5/thomas/Data/0503_Marseille/TOM03oct05
%         end
%     case 'fl'
%         catflag = true;
%         datatype = 2;
%         cd C:\Users\Thomas\WORK\IRMf\0503_Marseille\data\catflow\
% end

% get headers
[nfrms,xs,ys,hsize,nstims]=get_head(fname,'nframesperstim,framewidth,frameheight,lenheader,nstimuli');

% headers only
if nargin==2 && any(strfind(varargin{1},'header'))
    Y = struct('nfrms',nfrms,'xs',xs,'ys',ys,'hsize',hsize,'nstims',nstims);
    return
end

% Input
% frames
if nargin>1 && ~ischar(varargin{1})
    frames = varargin{1};
    nt = length(frames);
else
    frames = 1:nfrms;
end
nt = length(frames);
% sub-window
if nargin>2 && ~ischar(varargin{2})
    if nargin>3 && ~ischar(varargin{3})
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
Y = loadsum_single(fname,frames-1,0,datatype,xs,ys,hsize,indices);

% normalize
if nargin>1 && any(strfind(varargin{end},'norm'))
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