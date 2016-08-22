function Y = fast_normalizeseq(Yname,kernel1,kernel2)
% function Y = fast_normalizeseq(Y,kernel1[,kernel2])
% function Y = fast_normalizeseq('Y',kernel1[,kernel2])
%---
% new y(t) is calculated as kernel1*y(t) ./ kernel2*y(t-1), where * is the
% convolution
% default kernel1 is [1 1 1 1 1]
% default kernel2 is equal to kernel1

if nargin==0, help fast_normalize, end

% passage par r�f�rence �trange -> save memory !!
if ischar(Yname)
    Y = evalin('caller',Yname);
    evalin('caller',['clear ' Yname])
else
    Y = Yname;
end

% Input
if nargin<2
    kernel1 = [1 1 1 1 1];
else
    kernel1 = kernel1(:)';
end
if nargin<3
    kernel2 = kernel1;
else
    kernel2 = kernel2(:)';
end

% Normalize
[ni nj nt] = size(Y);
Y = reshape(Y,ni*nj,nt);
warning off MATLAB:divideByZero
disp('convolution1'), drawnow
if isscalar(kernel1)
    if kernel1~=1, error('scalar kernel should be equal to 1'), end
    Y1 = Y;
else
    Y1 = imfilter(Y,kernel1);
end
disp('convolution2'), drawnow
if isscalar(kernel2)
    if kernel2~=1, error('scalar kernel should be equal to 1'), end
    Y2 = Y(:,[2 end end]);
else
    kernel2bis = [kernel2 0 0]; % apply on t-1 instead of t
    Y2 = imfilter(Y,kernel2bis);
end
disp('divide')
Y = Y1./Y2;
warning on MATLAB:divideByZero
Y = reshape(Y,ni,nj,nt);
