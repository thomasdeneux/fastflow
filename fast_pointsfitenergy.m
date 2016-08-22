function e = fast_pointsfitenergy(iso,x,y,nj,ni)
% function e = fast_pointsfitenergy(iso,x,y,nj,ni)
%---
% iso is 1x3
% x and y must be 2xN

x(3,:)=1;
M = fast_register(iso,nj,ni);
xx = M*x;

diff = y-xx;
e = sum(diff(:).^2);