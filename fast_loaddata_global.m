function fast_loaddata_global(fname)
% function fast_loaddata_global(fname)
%---
% Load data from file fname into global variable Y. If Y already exists:
% - the class of Y (single or double) is preserved
% - if Y is already of the good size, the same memory space is used
%
% See also fast_loaddata, fast_loaddata_single

if nargin==0
    help fast_loaddata_global
    return
end

disp('assumed datype = 0')
datatype = 0;

global Y 

% get headers
[nfrms,xs,ys,hsize,nstims]=get_head(fname,'nframesperstim,framewidth,frameheight,lenheader,nstimuli');

% check size of Y
if ~isequal(size(Y),[xs ys nfrms])
    disp('new allocation for global array ''Y''')
    if isnumeric(Y), c = class(Y); else c = 'single'; end
    Y = []; %#ok<NASGU> (release memory space)
    Y = zeros([xs ys nfrms],c);
else
    disp('load data in existing global array ''Y''')
end

% load data
loadsum2_global(fname,datatype,xs,ys,nfrms,hsize);

