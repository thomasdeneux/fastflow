function Y = fast_load2bytes(fname,varargin)
% function Y = fast_load2bytes(fname[,Yname][,frames][,classname])
%---
% See also fast_save2bytes

fid = fopen(deblank(fname));
nj = fread(fid,1,'double');
ni = fread(fid,1,'double');
nt = fread(fid,1,'double');
fact = fread(fid,1,'double');

% passage par r�f�rence �trange -> save memory !!
Y = []; frames = 1:nt; classname = 'double';
for i=1:length(varargin)
    a = varargin{i};
    if isnumeric(a)
        frames = a;
    elseif regexp(a,'char|int|single|float|double')
        classname = a;
    else
        Y = evalin('caller',a);
        evalin('caller',['clear ' a])
        if any(size(Y)~=[nj ni length(frames)]), error('dimension mismatch'), end
    end
end
if isempty(Y)
    Y = zeros(nj,ni,length(frames),classname);
else
    classname = class(Y);
end

readflag = ['integer*2=>' classname];
for i=frames
    Yt = fread(fid,nj*ni,readflag);
    Y(:,:,i) = reshape(Yt/fact,nj,ni); %#ok<AGROW>
end

fclose(fid);