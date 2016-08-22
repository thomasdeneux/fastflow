function Z = fast_applymask(Y,mask)
% function Z = fast_applymask(Y,mask)

if nargin==0, help fast_applymask, end

[nj ni nt] = size(Y);
Z = reshape(Y,nj*ni,nt);
Z(~mask,:)=0;
Z = reshape(Z,nj,ni,nt);