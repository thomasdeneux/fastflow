function Y = fast_normalize(Yname)
% function Y = fast_normalize(Y)
% function Y = fast_normalize('Y')
%---
% achieves vessel registration (estimation of movement based on CS)
% and normalization (divide by mean)

if nargin==0, help fast_normalize, end

% passage par r�f�rence �trange -> save memory !!
if ischar(Yname)
    Y = evalin('caller',Yname);
    evalin('caller',['clear ' Yname])
else
    Y = Yname;
end

fprintf('dividing:           ')
warning off MATLAB:divideByZero
[nj ni nt] = size(Y);
Y = reshape(Y,nj*ni,nt);
M = mean(Y,2);
for i=1:nj*ni
    if ~mod(i,nj), fprintf('\b\b\b\b\b\b\b\b\b\b%4i/%4i\n',i/nj,ni), end
    Y(i,:) = Y(i,:)/M(i);
end
Y = reshape(Y,nj,ni,nt);
warning on MATLAB:divideByZero

