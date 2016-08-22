function mik_saverawfloat(Y,fname)
% function mik_saverawfloat(Y,fname)
%---
% See also mik_loadrawfloat

if nargin<1
    fname = fn_savefile('*.rawfloat');
end

[p b] = fileparts(fname);
fname = fullfile(p,[b '.rawfloat']);

fid = fopen(fname,'w');
fwrite(fid,size(Y),'uint32');
fwrite(fid,Y,'float32');
fclose(fid);
