function profile = fast_vesselwidth(u,p0)
% function profile = fast_vesselwidth(x[,p0])
%---
% Fit of the vessel section by a Gaussian + a linear trend
%
% Input:
% - x   np x nt array, each column is a vessel section at a different frame
% - p0  5 elements vector, initialization
%
% Output:
% - profile     5 x nt array, Gaussian parameters (width, center,
%               amplitude, baseline, trend)

u = double(1-u/nmean(u(:)));
[np nt] = size(u);
profile = zeros(5,nt);
% estimate: sigma, center, amplitude, baseline, trend
xi   =     [np/6   np/2    .1         prctil(u(:,1),10)  0];
xmin =     [1      np/4    0          -1        -1];
xmax =     [np/4   3*np/4  2          1         1];
if nargin==2, xi=double(p0); end
% options
opt = optimset('display','none','GradObj','on');
% binning
nbin = min(100,nt);
binsize = nt/nbin;
% loop
fn_progress('vessel width',nbin)
for i=[1:nbin nbin-1:-1:1]
    fn_progress(i)
    
    binidx = 1+round((i-1)*binsize):round(i*binsize);
    ui = mean(u(:,binidx),2);
    
    xi = fmincon(@(x)gaussianfit_(x,ui),xi,[],[],[],[],xmin,xmax,[],opt);
    profile(:,binidx) = repmat(xi(:),[1 length(binidx)]);
    
    
    % display
    if i==1
        hf = 579;
        if ~ishandle(hf)
            figure(hf)
            set(hf,'name','vessel width','numbertitle','off')
            axes
        end
        ha = findobj(hf,'type','axes');
        xcell = num2cell(xi);
        [sigm cent amplitude baseline trend] = deal(xcell{:});
        t = (1:np)';
        background = baseline + trend*(t/np-1/2);
        gaussian0  = exp(-(t-cent).^2/(2*sigm^2));
        gaussian   = amplitude * gaussian0;
        shape = background + gaussian;
        plot([ui shape],'parent',ha), pause(.005)
    end

end


%---
function [e de] = gaussianfit_(x,ui)

xcell = num2cell(x);
[sigma center amplitude baseline trend] = deal(xcell{:});

ok = ~isnan(ui);
np = length(ui);
t = (1:np)';

background = baseline + trend*(t/np-1/2);
gaussian0  = exp(-(t-center).^2/(2*sigma^2));
gaussian   = amplitude * gaussian0;
shape = background + gaussian;
error = ui - shape;

e = sum(error(ok).^2);

% figure(99)
% set(99,'tag','no-fn_imvalue')
% plot([ui shape]), title(num2str(e)), pause(.001)

if nargout==2
    dshape = [ ...
        (t-center).^2/sigma^3 .* gaussian ...
        (t-center)/sigma^2 .* gaussian ...
        gaussian0 ...
        ones(np,1) ...
        (t/np-1/2) ...
        ];
    de = -2*(error(ok)'*dshape(ok,:));
end
