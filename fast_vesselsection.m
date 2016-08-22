function varargout = fast_vesselsection(varargin)
% function par = fast_vesselsection('par') 
% function x = fast_vesselsection(Y,E,par)

switch nargin
    case 0
        help fast_vesselsection
        return
    case 1
        if ~strcmp(varargin{1},'par'), error argument, end
        varargout = {defaultparameters};
    case 3
        [Y E par] = deal(varargin{:});
        varargout = {vesselsection(Y,E,par)};
    otherwise
        error argument
end

%---
function par = defaultparameters

par.pos = .5;
par.dx = .8;
par.extentmode = 'relative';
par.extent = 4; % size section, as a ratio from vessel width
par.thicknessmode = 'pixel';
par.thickness = 4;

%---
function x = vesselsection(Y,E,par)

% section axis
idx = round(1+par.pos*(E.np-1));
c = E.snake(:,idx);
u = E.tubea(:,idx)-c;
un = u/norm(u);
switch par.extentmode
    case 'pixel'
        halfextent = par.extent / 2;
    case 'relative'
        halfextent = par.extent * norm(u);
    otherwise
        error('unknown extent mode ''%s''',par.extentmode)
end
segment = [c-halfextent*un c+halfextent*un];
sectionaxis = interp1([-halfextent halfextent],segment', ...
    -halfextent:par.dx:halfextent)';

% section thickness
switch par.thicknessmode
    case 'pixel'
        halfthickness = par.thickness / 2;
    case 'relative'
        halfthickness = par.thickness * norm(u);
    otherwise
        error('unknown thickness mode ''%s''',par.thicknessmode)
end
np = size(sectionaxis,2);
sectionhalfthickness = halfthickness * ones(1,np);

% interpolation
F = fast_edge(E.S,sectionaxis,sectionhalfthickness,par.dx);
x = interptube(F,Y);

