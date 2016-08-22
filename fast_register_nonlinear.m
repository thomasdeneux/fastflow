function varargout = fast_register_nonlinear(varargin)
% function [x frm2] = fast_register_nonlinear(ref,frm[,x0])
% function frm2 = fast_register_nonlinear(frm,x|xplus)
% function xplus = fast_register_nonlinear(frm,x,'parameters')
%---
% x contains a minimum-size description of the registration
% xplus contains a larger description which will improve resampling speed
% (usefull if the resampling is repeated many times with the same
% parameters)

% constants
gridsize = 20;
maxd = 10;
sigma = 1;
bin = ceil(gridsize/2);

% input
if isequal(size(varargin{2}),size(varargin{1}))
    doestimate   = true;
    doresample   = (nargout>=2);
    doparameters = doresample;
    [ref frm] = deal(varargin{1:2});
    if nargin>=3
        x = varargin{3};
    else
        x = [];
    end
else
    doestimate   = false;
    [frm tmp] = deal(varargin{1:2});
    doparameters = ~isstruct(tmp);
    if doparameters
        x = tmp;
    else
        xplus = tmp;
    end
    doresample = (nargin<3);
    if ~doresample && ~any(findstr(varargin{3},'parameters')), error('wrong flag'), end
end
dosparse = doparameters && ~doresample;
    
% general precomp
if doestimate || doparameters
    [ni nj nfr] = size(frm); %#ok<NASGU>
    [ii jj] = ndgrid(1:ni,1:nj);
    ngridx = max(ceil(ni/gridsize),2);
    ngridy = max(ceil(nj/gridsize),2);
    xgrid = linspace(1,ni,ngridx);
    ygrid = linspace(1,nj,ngridy);
    [xgrid ygrid] = ndgrid(xgrid,ygrid);
    if isempty(x)
        x = zeros(ngridx,ngridy,2);
    elseif numel(x) == ngridx*ngridy*2;
        x = reshape(x,ngridx,ngridy,2);
    else
        error('size mismatch for registration parameters')
    end
    nx = numel(x);
end

% estimation
if doestimate
    disp('nonlinear register'), tic
    if any(size(ref)~=[ni nj])
        error('size mismatch between frame and reference')
    end
    % precomp for energy
    ref1 = filt2(ref,sigma);
    frm1 = filt2(frm,sigma);
    subi = 2+maxd:ni-maxd-1;
    subj = 2+maxd:nj-maxd-1;
    [nibin njbin] = size(fn_bin(frm(subi,subj),[bin bin]));
    % prepare display
    hf = 328;
    if ~ishandle(hf)
        figure(hf)
        set(hf,'numbertitle','off','name','nonlinear register')
        ha = axes('parent',hf); 
        uicontrol('style','togglebutton','string','show figure', ...
            'parent',hf,'position',[10 10 80 20],'value',1);
    else
        ha = findobj(hf,'type','axes');
        hu = findobj(hf,'type','uicontrol');
        if get(hu,'value'), figure(hf), end
    end
    him = imagesc((frm1-ref1)','parent',ha);
    set(ha,'climmode','manual')
    % precomp for derivative
    dwarpidx = zeros(ni,nj,nx);
    dwarpjdx = zeros(ni,nj,nx);
    for k=1:nx
        xk = zeros(ngridx,ngridy,2); xk(k) = 1;
        dwarpidx(:,:,k) = interpn(xgrid,ygrid,xk(:,:,1),ii,jj,'spline');
        dwarpjdx(:,:,k) = interpn(xgrid,ygrid,xk(:,:,2),ii,jj,'spline');
    end
    [dfrmdj dfrmdi] = gradient(frm1);
    % estimate
    opt = optimset('jacobian','on','hessian','off', ...
        'Display','none','TolX',.02, ...
        'MaxIter',500,'MaxFunEvals',5000);
    maxd = ones(ngridx,ngridy,2)*maxd;
    x = lsqnonlin(@energy,x,-maxd,maxd,opt);
    fprintf('\b (%.1fs)\n',toc)
end

% parameters
if doparameters
    xplus = makegrid(x,dosparse);
end

% resample
if doresample
    frm2 = resample(frm,xplus);
end

% output
if doestimate
    if doresample
        varargout = {x frm2};
    else
        varargout = {x};
    end
else
    if doresample
        varargout = {frm2};
    else
        varargout = {xplus};
    end
end

% sub-functions
    function xplus = makegrid(x,dosparse)
        x = reshape(x,ngridx,ngridy,2);
        warpi = interpn(xgrid,ygrid,x(:,:,1),ii,jj,'spline');
        warpj = interpn(xgrid,ygrid,x(:,:,2),ii,jj,'spline');
        ii2 = warpi+ii;         jj2 = warpj+jj;
        if dosparse
            % compute a sparse array: more time to assemble it, but faster
            % resampling thereafter (useful if same resampling is applied
            % many times)
            okinterp = (ii2>1 & ii2<ni & jj2>1 & jj2<nj);
            okidx    = find(okinterp);
            ii3 = ii2(okinterp);    jj3 = jj2(okinterp);
            ii4 = floor(ii3);       jj4 = floor(jj3);
            u = ii3 - ii4;          v = jj3 - jj4;
            a = sparse( ...
                [okidx          okidx            okidx      okidx], ...
                [ii4+ni*(jj4-1) ii4+1+ni*(jj4-1) ii4+ni*jj4 ii4+1+ni*jj4], ...
                [(1-u).*(1-v)   u.*(1-v)         (1-u).*v   u.*v], ...
                ni*nj,ni*nj);
            xplus = struct('a',a,'ni',ni,'nj',nj);
        else
            % interpolation parameters only
            xplus = struct('ii2',ii2,'jj2',jj2);
        end
    end

    function frm2 = resample(frm,xplus)
        if isfield(xplus,'a')
            % sparse array multiplication
            frm2 = reshape(xplus.a*frm(:),xplus.ni,xplus.nj);
        else
            % interpolation - slower
            frm2 = interpn(frm,xplus.ii2,xplus.jj2); 
        end
    end

    function [e dedx frm2] = energy(x)
        % resample
        x = reshape(x,ngridx,ngridy,2);
        warpi = interpn(xgrid,ygrid,x(:,:,1),ii,jj,'spline');
        warpj = interpn(xgrid,ygrid,x(:,:,2),ii,jj,'spline');
        frm2 = interpn(frm1,ii+warpi,jj+warpj);
        frm2(isnan(frm2)) = mean(frm1(:));
        % energy: bin the square difference to avoid too large arrays
        d = frm2-ref1; 
        d2 = d.^2;
        e = fn_bin(d2(subi,subj),[bin bin]);
        % display
        set(him,'cdata',d'), drawnow
        % derivative        
        if nargout>=2
            dfrm2dwarpi = interpn(dfrmdi,ii+warpi,jj+warpj);
            dfrm2dwarpj = interpn(dfrmdj,ii+warpi,jj+warpj);
            dfrm2dwarpi(isnan(dfrm2dwarpi)) = 0;
            dfrm2dwarpj(isnan(dfrm2dwarpj)) = 0;
            dfrm2dx = repmat(dfrm2dwarpi,[1 1 nx]) .* dwarpidx ...
                + repmat(dfrm2dwarpj,[1 1 nx]) .* dwarpjdx;   
            dd2dx = 2 * repmat(d,[1 1 nx]) .* dfrm2dx;
            % energy derivative: bin too large arrays 
            dedx = fn_bin(dd2dx(subi,subj,:),[bin bin 1]);
            dedx = reshape(dedx,nibin*njbin,nx);
        end
    end

end