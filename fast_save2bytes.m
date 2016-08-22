function fast_save2bytes(Y,fname)
% function fast_save2bytes(Y,fname)
%
% See also fast_load2bytes

if nargin==0, help fast_save2bytes, return, end

[nj ni nt] = size(Y);

%disp('compute max')
m = max(max(Y(:)),-min(Y(:)));
fact = 2^15/m;

fid = fopen(fname,'w');
fwrite(fid,[nj ni nt fact],'double');

%fn_progress('saving:',nt)
for i=1:nt
    %fn_progress(i)
    Yt = Y(:,:,i);
    Yt = fix(Yt*fact);
    fwrite(fid,Yt,'integer*2');
end

fclose(fid);
