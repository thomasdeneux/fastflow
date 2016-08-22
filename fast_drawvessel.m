function fast_drawvessel(e,CS)
% function fast_drawvessel(e[,CS])

if nargin==0, help fast_drawvessel, end

if nargin>1
    imagesc(CS)
end
for i=1:numel(e)
    x = e(i).points2;
    line(x(:,1),x(:,2))
end