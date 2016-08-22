function Y = fast_loadrawfloat(fname)
% function fast_loadrawfloat(Y,fname)
%---
% See also fast_saverawfloat

if nargin<1
    fname = fn_getfile('*.rawfloat');
end

fid = fopen(fname,'r');
s = fread(fid,3,'uint32')';
Y = reshape(fread(fid,prod(s),'uint32'),s);
fclose(fid);
