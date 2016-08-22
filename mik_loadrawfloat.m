function Y = mik_loadrawfloat(fname)
% function mik_loadrawfloat(Y,fname)
%---
% See also mik_saverawfloat

if nargin<1
    fname = fn_getfile('*.rawfloat');
end
[p b] = fileparts(fname);
fname = fullfile(p,[b '.rawfloat']);

fid = fopen(fname,'r');
s = fread(fid,3,'uint32')';
Y = reshape(fread(fid,prod(s),'float32'),s);
fclose(fid);
