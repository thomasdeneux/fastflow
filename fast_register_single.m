function out = fast_register_single(varargin)
% same as fast_register, but heavy computations and are done in class
% single, and CSinv output is of class single
%---
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
% - CSinv   CS(iso^-1(x,y))  -  class single 
% - M       associated matrix (3*4)

if nargin==0, help fast_register_single, return, end 

% 1) function M = fast_register(iso,nj,ni)
if length(varargin{2})==1
    if any(size(varargin{1})==1)
        [iso nj ni] = deal(varargin{:});
        out = makematrice(iso,nj,ni);
        return
% 2) function iso = fast_register(M,nj,ni)
    else
        [M nj ni] = deal(varargin{:});
        out = findiso(M,nj,ni);
        return   
    end
end

% 3) function CSinv = fast_register(CS,iso)
if length(varargin{2})==3
    [CS iso] = deal(varargin{:});
    %out = invisometrie(CS,iso);
    out = fast_regenergy_single(iso,single(CS));
    return
end

% 4) function isometrie = fast_register(CS1,CS2[,iso0[,indices]])
[CS1 CS2] = deal(varargin{1:2});
CS1 = single(CS1);
CS2 = single(CS2);
if nargin>2, iso0 = varargin{3}; iso0=iso0(:)'; else iso0 = [0 0 0]; end
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
   'Display','none','TolX',1e-6,'TolFun',1e-6,'MaxIter',maxiter,'MaxFunEvals',maxiter*100);

allowrotations = true;
% when calculating derivatives, differences are too small once double are
% converted to single! the 'relative' difference (i.e. deltax/x) can be
% artificially increased by making the value to be estimated in lsqnonlin
% orders of magnitude smaller
fact = 1e5; 
if allowrotations
    %out = lsqnonlin(@energy,iso0,[],[],opt,CS1,CS2,indices);
    out = fact * lsqnonlin(@test,iso0/fact,[],[],opt,fact,CS1,CS2,indices);
else
    % on n'estime que les 2 param�tres de translation ; la fonction
    % fast_regenergy.dll rep�re toute seule que s'il y a 3 param�tres c'est
    % [rotation transx transy], et s'il y en a seulement 2, c'est [transx
    % transy]
    
    %out = lsqnonlin(@energy,iso0(2:3),[],[],opt,CS1,CS2,indices);
    out = fact * lsqnonlin(@test,iso0(2:3)/fact,[],[],opt,CS1,CS2,indices);
    out = [0 out];
end

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
function e = test(x,fact,   varargin)

e = fast_regenergy_single(x*fact,varargin{:});
