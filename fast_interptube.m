function y = fast_interptube(Y,snake,halfd,tubea,tubeb)
% function y = fast_interptube(Y,snake,halfd,tubea,tubeb)
% function y = fast_interptube(Y,edge)

if nargin==2
    e = snake;
    snake = e.snake;
    halfd = e.halfd;
    tubea = e.tubea;
    tubeb = e.tubeb;
end
np = size(snake,2);

% maximum size for the filter: idx of pixel in the
% center
icenter = ceil(max(halfd));
nf = 2*icenter;
% coordinates of points in the grid
[xx yy] = ndgrid((1:nf)-icenter,(1:nf)-icenter,ones(1,np));
% coordinates of points in image - dimensions are
% (gridi,gridj,kpoint,xy)
c = repmat(shiftdim(floor(snake'),-2),nf,nf);
x = cat(4,xx,yy) + c;
% distance to point A, point B and line
a = repmat(shiftdim(tubea',-2),nf,nf);
b = repmat(shiftdim(tubeb',-2),nf,nf);
ax = x-a;
bx = x-b;
da = sqrt(sum(ax.^2,4));
db = sqrt(sum(bx.^2,4));
ab = b-a;
dab = sqrt(sum(ab.^2,4));
u = ab ./ repmat(dab,[1 1 1 2]);
cross = u(:,:,:,1).*ax(:,:,:,2) - u(:,:,:,2).*ax(:,:,:,1);
dd = abs(cross);
% distance to segment
dot = sum(u.*ax,4);
left = (dot<0);
right = (dot>dab);
dist = left.*da + (~left & ~right).*dd + right.*db;
% filter
H = max(1-dist,0);
H = H ./ repmat(sum(sum(H)),nf,nf);

% parameters for interpolation
p = repmat(shiftdim(1:np,-1),nf,nf);
f = logical(H);
A.p = p(f);
A.i = c(:,:,:,1)+xx; A.i = A.i(f);
A.j = c(:,:,:,2)+yy; A.j = A.j(f);
A.z = H(f);
% build sparse interpolation matrix
np = size(snake,2);
[ni nj nt] = size(Y);
a = sparse(A.p,double(A.i+ni*(A.j-1)),double(A.z),np,ni*nj);
% interpolation
Y = reshape(Y,ni*nj,nt);
y = a*double(Y);

