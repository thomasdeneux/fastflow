function varargout = fast_edges2image(varargin)
% function s = fast_edges2image(S,'matrix')
% function y = fast_edges2image(x,s,background)
%---
% Interpolation from vessel content to the full image
%
% Input/Output:
% - S           fastflow object
% - 'matrix'    flag to get the sparse interpolation array
% - s           parameters for interpolation
% - x           cell array of 2D (or 3D for colors) arrays, x{k} must be of
%               size np(k) x n 
% - background  background image (color image if x is 3D)
% - y           output movie

if nargin==0, help fast_edges2image, return, end

if ischar(varargin{2}) && strcmp(varargin{2},'matrix')
    varargout = {getsparsearray(varargin{1})};
else
    varargout = {interpolate(varargin{:})};
end

%---
function s = getsparsearray(S)

edges  = S.edges([S.edges.active]);
a = sparse(S.nx*S.ny,0);
for i=1:numel(edges)
    e = edges(i);
    x=double([e.tubea e.tubeb(:,end:-1:1)]);
    mask=poly2mask(x(2,:),x(1,:),S.nx,S.ny);
    [ii jj] = find(mask); % column vectors
    dd = fn_add(ii,-e.snake(1,:)).^2 + fn_add(jj,-e.snake(2,:)).^2;
    [dum kk] = min(dd,[],2); % indices in vessels corresponding to minimal distance
    ai = sparse(ii+S.nx*(jj-1),kk,1,S.nx*S.ny,e.np);    
    a = [a ai]; %#ok<AGROW>
end

indices = any(a,2);
a = a(indices,:);
np = sum(indices);
a = spdiags(1./sum(a,2),0,np,np)*a;
s = struct('a',a,'nx',S.nx,'ny',S.ny,'indices',indices);

%---
function y = interpolate(x,s,background)

if isscalar(background), background = background*ones(s.nx,s.ny); end

if iscell(x), x = cat(1,x{:}); end
x = double(x);
nfr = size(x,2);
nc = size(x,3); if ~ismember(nc,[1 3]), error dimension, end
y = repmat(fn_reshapepermute(background,{[1 2] 4 3}),[1 nfr]);
for k=1:nc, y(s.indices,:,k) = s.a*x(:,:,k); end
y = reshape(y,[s.nx s.ny nfr nc]);



