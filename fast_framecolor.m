function frs = fast_framecolor(y,yclip,x,a,xclip,xcmap,vesseltype)
% function frs = fast_framecolor(y,yclip,x,a,xclip,xcmap,vesseltype)
%---
% Input:
% - y       3D (nx x ny x nt) array - movie
% - yclip   clipping for movie
% - x       2D (np x nt) array or cell array - movie of vessels
% - a       sparse array - interpolation of the vessels
% - xclip   clipping for vessels
% - xcmap   colormap for vessels
% - vesseltype  mask with 0 outside vessels, 1 for arteries, 2 for veins
% Output:
% - frs     4D (nx x ny x 3 x nt) array - movie with vessels inside

% Data
if iscell(x)
    x = cat(1,x{:});
end
[nx ny nt] = size(y);
y = fn_clip(y,yclip);
if nt==1
    nt = size(x,2);
    y = repmat(y,[1 1 nt]);
end
np = size(x,1); 
if size(x,2)~=nt || size(a,1)~=nx*ny || size(a,2)~=np
    error('size mismatch')
end

% Clipping + colormap
x = fn_clip(x,xclip,xcmap); 
x = reshape(permute(x,[1 3 2]),[np 3*nt]);

% Vessels type
y = reshape(y,[nx*ny 1 nt]);
if nargin>=7
    maskartery = bwmorph(vesseltype==1,'dilate',0);
    maskvein   = bwmorph(vesseltype==2,'dilate',0) & ~maskartery;
    mask = .9*maskartery + .5*maskvein;
    bmask = logical(mask);
    y(bmask,:,:) = repmat(mask(bmask),[1 1 nt]);
end


% Movie
frs = repmat(y,[1 3]);
vesselmask = any(a,2);
frs(vesselmask,:) = a(vesselmask,:)*x;
frs = reshape(frs,[nx ny 3 nt]);

