function varargout = fast_fakedata(varargin)
% function [y pos] = fast_fakedata(par)
% function par = fast_fakedata('par')
% function par = fast_fakedata('userpar'[,par])
%---
% Create simulated data of RBCs moving in a vessel.
% 
% Input:
% - par     structure with parameters:
%           .step       speed of the RBCs
%           .density    probability of a RBC at a given position
%           .diam       RBCs diameter
%           .jitter     amplitude of random fluctuation in speed (note that
%                       jitter is NOT scaled by step)
%           .shotnoise  amplitude of additional shot noise
%           .change     relative change in speed for the middle part of the
%                       time segment
%           .svx        spatial smoothness of the speed
%           .svt        temporal smoothness of the speed
%           .nx         number of points of the blood vessel
%           .nt         number of time instants
%
% Output:
% - y       nx*nt array: lines picture
% - pos     nRBC*nt array: positions of the RBCs accross time
% 
% Ivo, 23 June 2006nargin==1
% Thomas, 26 October 2010

% Input
if ischar(varargin{1})
    switch(varargin{1})
        case 'par'
            par = defaultparameters;
        case 'userpar'
            par = defaultparameters;
            if nargin>=2
                par = fn_structmerge(par,varargin{2},'skip');
            end
            par = userparameters(par);
        otherwise
            error('wrong flag ''%s''',varargin{1})
    end
    varargout = {par};
else
    par = varargin{1};
    par = fn_structmerge(defaultparameters,par,'skip');
    nout = max(1,nargout);
    varargout = cell(1,nout);
    [varargout{:}] = createfake(par);
end

%---
function [par spec] = defaultparameters

par.step        = 1;
par.change      = .05;
par.density     = .1;
par.diam        = 3;
par.jitter      = 1;
par.svx         = .5;
par.svt         = 10;
par.shotnoise   = .1;
par.nx          = 100;
par.nt          = 300;
par.nexp        = 512;

spec.step       = 'stepper 1 0 20 .5 %.3g';
spec.change     = 'stepper 1 0 .5 .01 %.3g';
spec.density    = 'stepper 1 .01 .5 .05 %.3g';
spec.diam       = 'stepper 1 1 30 .5 %.3g';
spec.jitter     = 'stepper 1 0 100 1 %.2g';
spec.svx        = 'stepper 1 0 1 .1 %.3g';
spec.svt        = 'stepper 1 0 100 1 %i';
spec.shotnoise  = 'stepper 1 0 10 .1 %.2g';
spec.nx         = 'stepper 1 10 200 10 %i';
spec.nt         = 'stepper 1 10 2000 100 %i';
spec.nexp       = 'stepper 1 1 Inf 2 %i';

%---
function par = userparameters(par)

[dum spec] = defaultparameters;
hf = figure(942); clf(hf)
set(hf,'numbertitle','off','name','fake parameters', ...
    'defaultuicontrolfontsize',8, ...
    'position',[300 200 800 420])
hp = uipanel('parent',hf,'position',[.01 .1 .2 .8],'units','pixel');
ha = axes('parent',hf,'position',[.3 .1 .68 .87]);
ok = uicontrol('parent',hf,'style','togglebutton','string','ok', ...
    'units','normalized','position',[.93 .01 .06 .04]);
X = fn_control(par,spec,@(s)showfake(s,ha),hp);
showfake(par,ha)
waitfor(ok,'value',1)
if isvalid(X), par = X.s; else par=[]; end
if ishandle(hf), close(hf), end

%---
function showfake(par,ha)

y = createfake(par);
imagesc(y,'parent',ha)

%---
function [y pos] = createfake(par)

% check
if par.step<0 || par.change<0, error('step and change must be non-negative'), end

% parameters
F0      = 100; % baseline signal
RBCamp  = -1;  % amplitude (contrast) of RBCs 
sigma   = par.diam/2;  % lateral fall-off of RBCs (i.e. RBCs are modelled as a gaussian)


% creation of baseline signal and shot noise
y = F0 + par.shotnoise * randn(par.nx,par.nt);

% initial distribution of RBC positions (a vector with values between
% -nresv-nresjit and +nx+nresjit)
nreservev   = ceil(par.nt*par.step*(1+par.change));
nreservejit = 2*ceil(sqrt(par.nt)*par.jitter);
a    = rand(par.nx+nreservev+2*nreservejit,1);
pos0 = find(a<par.density) - (nreservev+nreservejit); % column vector
nRBC = length(pos0);

% velocity of the RBC
% (regular motion)
begchange = floor(par.nt/3);
endchange = floor(2*par.nt/3);
v = zeros(nRBC,par.nt);
v(:,[1:begchange endchange+1:par.nt]) = par.step;
v(:,begchange+1:endchange)        = par.step*(1+par.change);
% (jitter)
jitcommon = randn(1,par.nt) * (par.jitter*sqrt(par.svx));       % part of jitter common to all RBCs
jitindept = randn(nRBC,par.nt) * (par.jitter*sqrt(1-par.svx));  % part of jitter independent for each RBC
jit = fn_add(jitcommon,jitindept);
jit = filtx(jit,par.svt) * max(1,sqrt(par.svt)/3); % temporal smoothing and heuristic correction to preserve std
v = v+jit;

% position of RBCs (attention, need to keep RBCs sorted)
pos = fn_add(pos0,cumsum(v,2));

% 'print' each RBC in the line image
pos1     = permute(pos,[3 2 1]); % 1 x nt x nRBC array
dist     = fn_add((1:par.nx)',-pos1); % nx x nt x nRBC array - distance of all RBCs at all times from all pixels
contrast = RBCamp * exp(-(dist/sigma).^2/2);
y = y + sum(contrast,3);

