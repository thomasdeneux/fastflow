function fast_load2bytes_global(fname)
% function fast_load2bytes_global(fname)
%---
% See also fast_save2bytes

fid = fopen(fname);
nj = fread(fid,1,'double');
ni = fread(fid,1,'double');
nt = fread(fid,1,'double');
fact = fread(fid,1,'double');
frames = 1:nt;

global Y 

% check size of Y
if isempty(Y) || any(size(Y)~=[nj ni nt])
    disp('new allocation for global array ''Y''')
    c = class(Y);
    Y = []; %#ok<NASGU> (release memory space)
    Y = zeros([nj ni nt],c);
else
    disp('load data in existing global array ''Y''')
end

% read
readingmode = ['integer*2=>' class(Y)];
for i=frames
    Yt = fread(fid,nj*ni,readingmode);
    Y(:,:,i) = reshape(Yt/fact,nj,ni);
end

fclose(fid);