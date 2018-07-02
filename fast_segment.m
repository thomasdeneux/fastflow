function [out evol] = fast_segment(varargin)
% function edge = fast_segment([mask,]poly|CS[,ha][,par])
% function par = fast_segment('par')
% function paramstr = fast_segment('parameters')
%---
% Input
% - mask    image of "energy field" (should be positive inside the vessels,
%           and negative outside of the vessels, typically use mask =
%           fast_structure(CS), where CS is the reference vasculature
%           image)
% - image   raw data image
% - poly    2xn array - defines the user initialization of the vessel
% - CS      image - if specified (alternatively to 'poly'), the user is
%           prompted to specify this initialization with the mouse
% - ha      vector of axes handles where to show the convergence
% - par     structure with estimation parameters
%
% Output
% - edge    structure with fields 'snake' (skeleton of the vessel),
%           'halfd' (half diameter), 'tubea' and 'tubeb' (sides)
% - paramstr  returns typical parameter values in the form of a character
%           array which can be evaluated

if nargin==1 && ischar(varargin{1})
    switch varargin{1}
        case 'parameters'
            out = defaultpar_string;
        case 'par'
            out = defaultpar;
        otherwise
            error('unknonwn flag "%s"', varargin{1})
    end
    return
else
    if nargin<2 || (ischar(varargin{2}) || isstruct(varargin{2}))
        % no energy mask provided, compute it
        x = double(varargin{1}); % image
        mask = fast_structure(x);
        karg = 2;
    else
        mask = varargin{1};
        x = varargin{2};        % image or polygon initialization
        karg = 3;
    end
end

mask = mask/nstd(mask(:)); % normalize energy field by STD
par = defaultpar;
hlist = [];
for i=karg:length(varargin)
    a = varargin{i};
    if iscell(a) || isstruct(a)
        par = a;
    elseif ischar(a) 
        par = cellstr(a);
    else
        hlist = a(:)'; % list of axes handles
    end
end

if size(x,1)==2
    poly = x;
elseif size(x,2)==2;
    poly = x';
else
    CS = x;
    hf = fn_figure('FAST_SEGMENT');
    figure(hf)
    fn_imvalue image
    subplot(121), imagesc(CS')
    subplot(122), imagesc(mask')
    h = msgbox({'Select region, press OK,','then select vessel initialization'});
    waitfor(h)
    poly = fn_mouse('poly');
    if size(poly,1)~=2, error programming, end
    hlist(end+1) = gca;
end

% parameters
if iscell(par)
    for i = 1:length(par)
        eval(par{i});
    end
else
    % structure
    F = fieldnames(par);
    for k=1:length(F)
        f = F{k}; val = par.(f); %#ok<NASGU>
        eval([f '=val;']);
    end
end

% ENERGY TO MAXIMIZE
% E = sum_insidex(mask)
% dE/dx = sum_x(mask*dxtowardsoutside)

% trick the energy field
f = (mask<0);
mask(f) = mask(f)*WNEG;

% initialization
dpoly = sqrt(sum(diff(poly,1,2).^2));
poly(:,~dpoly) = []; % remove repeated points
dpoly(~dpoly) = [];
taille = [0 cumsum(dpoly)];
phi = 0:DX:taille(end);
np = length(phi);
snake = interp1(taille,poly',phi)';
halfd  = ones(1,np);

% display
hl = zeros(length(hlist),3);
for i = 1:length(hlist)
    axes(hlist(i))
    hl(i,1) = line('color','b','erasemode','xor','linestyle','--');
    hl(i,2) = line('color','b','erasemode','xor','linestyle','--');
    hl(i,3) = line('color','b','erasemode','xor');
end
setline(hl,snake,halfd)

% optimization
conv = false;
niter = 300;
if nargout>=2 % save the whole evolution?
    evol = zeros(3,np,niter+1); 
    evol(:,:,1) = [snake; halfd];
end
for i=1:niter
    dsnake = gradient(snake);
    normd  = sqrt(sum(dsnake.^2));
    
    % points equally spaced
    dnormd = gradient(normd);
    dsnake1 = dsnake .* repmat(dnormd,2,1);
    
    % snake smoothness
    d2snake = gradient(dsnake);
    % local component + regional component!
    dsnake2 = d2snake + smooth(d2snake,LSS);
    % ponderation by the width!
    dsnake2 = dsnake2.*repmat(halfd,2,1);
    
    % width smoothnes
    shalfd = smooth(halfd,LSD);
    dhalfd2 = shalfd - halfd;
    
    % energy
    [tubea tubeb] = gettube(snake,halfd);
    weighta = interpn(mask,tubea(1,:),tubea(2,:));
    weightb = interpn(mask,tubeb(1,:),tubeb(2,:));
    dtubea = weighta;
    dtubeb = weightb;
    dhalfd3 = (dtubea+dtubeb)/2;
    dnorm  = sqrt(sum(dsnake.^2));
    snaken = [dsnake(2,:)./dnorm; -dsnake(1,:)./dnorm];
    dsnake3 = snaken .* repmat((dtubea-dtubeb)/2,2,1);
    dsnake3 = smooth(dsnake3,LEN);
    
    dsnake = KEQ*dsnake1 + KSS*dsnake2 + KEN*dsnake3;    
    dhalfd = KSD*dhalfd2 + KEN*dhalfd3;
    
    % distance to anchors
    dsnake(:,1) = dsnake(:,1) + KEX*(poly(:,1)-snake(:,1));
    dsnake(:,end) = dsnake(:,end) + KEX*(poly(:,end)-snake(:,end));

    pause(.02)
    
    % update
    snake1 = snake + DT*dsnake;
    halfd1 = halfd + DT*dhalfd;
    halfd1 = min(max(halfd1,MINHALFD),MAXHALFD);
    
    % save the whole evolution?
    if nargout>=2
        evol(:,:,i+1) = [snake1; halfd1]; 
    end
    
    % divergence check
    if any(max(snake1,[],2)'+max(halfd1) > size(mask)) ...
            || any(min(snake1,[],2)'-max(halfd1) < [1 1])
        if nargout>=2, evol(:,:,i+2:end) = []; end
        break
    else
        snake = snake1;
        halfd = halfd1;
    end
    
    % display
    setline(hl,snake,halfd)
    
    % convergence detection
    normd = sqrt(sum(dsnake.^2));
    if max(normd)<CONVTHR
        conv = true; 
        if nargout>=2, evol(:,:,i+2:end) = []; end
        break
    end
end
if conv
    disp('convergence success')
else
    disp('convergence failed')
end

% resampling and narrowing tube
dsnake = sqrt(sum(diff(snake,1,2).^2));
taille = [0 cumsum(dsnake)];
phi = (0:DX2:taille(end))';
snake = interp1(taille,snake',phi)';
halfd = interp1(taille,halfd',phi)' * WFACT;

% output
out = struct('snake',snake,'halfd',halfd,'dx',DX2);

% clear display
delete(hl)

%---
function snaken = orthvect(snake)

dsnake = gradient(snake);
dnorm  = sqrt(sum(dsnake.^2));
snaken = [dsnake(2,:)./dnorm; -dsnake(1,:)./dnorm];

%---
function [tubea tubeb] = gettube(snake,halfd)

snaken = orthvect(snake);
halfd = repmat(halfd,2,1);
tubea = snake + snaken.*halfd;
tubeb = snake - snaken.*halfd;

%---
function x = smooth(x,L)

x = imfilter(x,fspecial('gaussian',[1 L*2],L),'replicate');

%---
function setline(hl,snake,halfd)

if isempty(hl), return, end

[tubea tubeb] = gettube(snake,halfd);

for i=1:size(hl,1)
    set(hl(i,1),'xdata',tubea(1,:),'ydata',tubea(2,:))
    set(hl(i,2),'xdata',tubeb(1,:),'ydata',tubeb(2,:))
    set(hl(i,3),'xdata',snake(1,:),'ydata',snake(2,:))
end

%---
function setarrows(x,dx)

delete(findobj(gca,'tag','setarrows'))
hold on
quiver(x(1,:),x(2,:),dx(1,:),dx(2,:),'tag','setarrows','erasemode','xor')
hold off

%---
function paramstr = defaultpar_string()

paramstr = {
    '% spacing between points'
    'DX = .5;'
    '% weights for evolution'
    'KEX = .02;    % extremities'
    'KEQ = .1/DX;  % equally spaced'
    'KSS = .2/DX;  % snake smoothing'
    'KSD = 1;      % vessel diameter smoothing'
    'KEN = 5e-3;   % energy factor'
    '% global factor for parameters above'
    'DT = 1;'
    '% smoothing parameters'
    'LSS = 5/DX;'
    'LSD = 20/DX;'
    'LEN = 5/DX;'
    '% discouragement for negative regions'
    'WNEG = 1;'
    '% constraints on vessel diameter'
    'MINHALFD = .5;'
    'MAXHALFD = 4;'
    '% convergence detection'
    'CONVTHR = 2e-2;'
    '% resampling'
    'DX2 = .8;'
    '% narrow tube'
    'WFACT = .7;'
    }; 

%---
function par = defaultpar

par = struct;
par.DX = .5;
par.KEX = .02;    % extremities
par.KEQ = .1/par.DX; % equally spaced
par.KSS = .2/par.DX; % smoothing
par.KSD = 1;     % smoothing
par.KEN = 5e-3;   % energy
par.DT = 1;
par.LSS = 5/par.DX;
par.LSD = 20/par.DX;
par.LEN = 5/par.DX;
par.WNEG = 1;
par.MINHALFD = .5;
par.MAXHALFD = 4;
par.CONVTHR = 2e-2;
par.DX2 = .8;
par.WFACT = .7;


