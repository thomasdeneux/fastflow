function out = fast_register(varargin)
% function isometrie = fast_register(CS1,CS2[,iso0[,indices]])
% function CSinv = fast_register(CS,iso)
% function M = fast_register(iso,nj,ni)
% function iso = fast_register(M,nj,ni)
%---
% Input:
% - CS1     reference image
% - CS2     target image
% - nj,ni   size of CS1
% - iso     isometrie
% Output:
% - isometrie CS1->CS2      [theta tx ty]
% - CSinv   CS(iso^-1(x,y))
% - M       associated matrix (3*4)

if nargin==0, help fast_register , return, end 

% function M = fast_register(iso,nj,ni)
if length(varargin{2})==1
    if any(size(varargin{1})==1)
        [iso nj ni] = deal(varargin{:});
        out = makematrice(iso,nj,ni);
        return
% function iso = fast_register(M,nj,ni)
    else
        [M nj ni] = deal(varargin{:});
        out = findiso(M,nj,ni);
        return   
    end
end

% function CSinv = fast_register(CS,iso)
if length(varargin{2})==3
    [CS iso] = deal(varargin{:});
    switch class(CS)
        case 'single'
            out = fast_regenergy_single(iso,CS);
        case 'double'
            out = fast_regenergy(iso,CS);
        otherwise
            error argument
    end
    return
end

% function isometrie = fast_register(CS1,CS2[,iso0[,indices]])
[CS1 CS2] = deal(varargin{1:2});
if nargin>2
    iso0 = varargin{3}; 
    iso0=iso0(:)'; 
    if isempty(iso0), iso0 = [0 0 0]; end
else
    iso0 = [0 0 0];
end
if nargin>3, 
    indices = varargin{4};
else
    [nj ni] = size(CS1);
    skipe = 20;
    [I,J] = meshgrid(1+skipe:ni-skipe,1+skipe:nj-skipe);
    indices = J(:) + nj*(I(:)-1);
end

% % function isometrie = fast_register(CS1,CS2[,iso0[,indices]])
% [CS1 CS2] = deal(varargin{1:2});
% if nargin>2, iso0 = varargin{3}; iso0=iso0(:)'; iso0(2:3)=iso0(2:3)/3; else iso0 = [0 0 0]; end
% if nargin>3, 
%     disp('don''t give indices, now a new method is tried using the whole image')
%     indices = varargin{4};
% end
% CS1 = CS1(2:3:end,2:3:end); CS2 = CS2(2:3:end,2:3:end);
% [nj ni] = size(CS1);
% skipe = 7;
% [I,J] = meshgrid(1+skipe:ni-skipe,1+skipe:nj-skipe);
% indices = J(:) + nj*(I(:)-1);
 
% maxiter = 500;
% opt = optimset('jacobian','off','largescale','off','LevenbergMarquardt','on', ...
%     'Display','off','TolX',1e-5,'TolFun',1e-5,'MaxIter',maxiter,'MaxFunEvals',maxiter*10);

maxiter = 500;
opt = optimset('jacobian','off', ...
   'Display','none','TolX',1e-6,'TolFun',1e-10,'MaxIter',maxiter,'MaxFunEvals',maxiter*100);

allowrotations = false;
switch class(CS1)
    case 'single'
        rescale = true; % single: need to rescale iso to avoid gradient zero
        fun = @regenergy1000;
        iso0 = iso0/1000;
    case 'double'
        rescale = false;
        fun = @fast_regenergy;
    otherwise
        error argument
end
if allowrotations
    %out = lsqnonlin(@energy,iso0,[],[],opt,CS1,CS2,indices);
    out = lsqnonlin(fun,iso0,[],[],opt,CS1,CS2,indices);
else
    % on n'estime que les 2 param???tres de translation ; la fonction
    % fast_regenergy.dll rep???re toute seule que s'il y a 3 param???tres c'est
    % [rotation transx transy], et s'il y en a seulement 2, c'est [transx
    % transy]
    
    %out = lsqnonlin(@energy,iso0(2:3),[],[],opt,CS1,CS2,indices);
    out = lsqnonlin(fun,iso0(2:3),[],[],opt,CS1,CS2,indices);
    out = [0 out];
end
if rescale, out = out*1000; end

% out(2:3) = out(2:3)*3;

%---
function M = makematrice(iso,nj,ni)

% VERSION C++ EGALEMENT

theta = iso(1);
shift = iso(2:3); shift=shift(:);
R = [cos(theta) -sin(theta) 
    sin(theta)  cos(theta)];
C = [ni/2 ; nj/2];
M = [R (shift + C - R*C)];

%---
function iso = findiso(M,nj,ni)

R = M(1:2,1:2);
theta = acos(R(1,1));
C = [ni/2 ; nj/2];
shift = M(1:2,3) - C + R*C;
iso = [theta shift'];

%---
function CSinv = invisometrie(CS,iso)

% VERSION C++ : fast_regenergy(iso,CS)

[nj ni] = size(CS);
M = makematrice(iso,nj,ni);

[I J] = meshgrid(1:ni,1:nj);
coords = [I(:)' ; J(:)'];
coords1 = [coords ; ones(1,prod(size(I)))];
coordsM = M * coords1;
I2 = reshape(coordsM(1,:),nj,ni);
J2 = reshape(coordsM(2,:),nj,ni);

CSinv = interp2(CS,I2,J2);

%---
%function [E, J] = energy(iso,CS1,CS2,indices)
function E = energy(iso,CS1,CS2,indices)

% VERSION C++ : fast_regenergy(iso,CS1,CS2,indices)

pause(0)
drawnow

[nj ni] = size(CS1);
M = makematrice(iso,nj,ni);

nindices = length(indices);
[I J] = meshgrid(1:ni,1:nj);
coords = [I(indices)' ; J(indices)'];
coords1 = [coords ; ones(1,nindices)];
coordsM = M * coords1;

CS2inv = interp2(CS2,coordsM(1,:),coordsM(2,:));
diff = CS2inv(:)-CS1(indices(:));
% % Output: one value
% E = diff(:)'*diff(:)/(nki*nkj)/1000;
% Output: the vector of all differences
E = diff(:)/100;
E(isnan(E)) = 1e8;

% gradient: NEEDS TO BE CHECKED AGAIN !!!
if nargout>1
    [dCS2dI dCS2dJ] = gradient(CS2);
    dCS2invdI = interp2(dCS2dI,I2,J2);
    dCS2invdJ = interp2(dCS2dJ,I2,J2);

    dRdtheta = [-sin(theta) -cos(theta) 
                cos(theta)  -sin(theta)];
    dMdtheta = [dRdtheta -dRdtheta*C];
    dcoordsMdtheta = dMdtheta * coords1;
    dI2dtheta = reshape(dcoordsMdtheta(1,:),nj-2*skip,ni-2*skip);
    dJ2dtheta = reshape(dcoordsMdtheta(2,:),nj-2*skip,ni-2*skip);
    dCS2dtheta = dCS2invdI.*dI2dtheta + dCS2invdJ.*dJ2dtheta;
    dEdtheta = 2*dCS2dtheta(:)'*diff(:);
    
    dCS2dshiftI = dCS2invdI;
    dEdshiftI = 2*dCS2dshiftI(:)'*diff(:);
    
    dCS2dshiftJ = dCS2invdJ;
   dEdshiftJ = 2*dCS2dshiftJ(:)'*diff(:);
    
    J = [dEdtheta dEdshiftI dEdshiftJ]/(ni*nj);
end

%---
function E = regenergy1000(iso,CS1,CS2,indices)

iso = iso*1000;
E = fast_regenergy_single(iso,CS1,CS2,indices);






