function varargout = fast_track(varargin)
% function par = fast_track('par') [get default parameters]
% function [J V] = fast_track(I,par[,J0,V0][,ha])
% function [J V] = fast_track('cut',I,par[,J0,V0])
% function J = fast_track('fakelines',V)
%---
% note that here, J0, V0, J and V have the same size as I

if ischar(varargin{1})
    % Special actions
    switch varargin{1}
        case 'cut'
            % Estimate smaller blocks separately
            BLOCKSIZE = 150;
            SIDE = 15;
            [I par] = deal(varargin{2:3});
            par = parameters(I,par);
            if nargin<5
                [J0 V0] = init(I,par);
            else
                [J0 V0] = deal(varargin{4:5});
            end
            nblocks = ceil(par.nt/BLOCKSIZE);
            blocksize = (par.nt-1)/nblocks;
            J = zeros(par.nx,par.nt);
            V = zeros(par.nx,par.nt);
            for k=1:nblocks
                idx0 = 1+(round((k-1)*blocksize):round(k*blocksize));
                idx = idx0;
                subidx = 1:length(idx);
                if k>1, idx = [idx(1)+(-SIDE:-1) idx]; subidx = subidx+SIDE; end %#ok<AGROW>
                if k<nblocks, idx = [idx idx(end)+(1:SIDE)]; end %#ok<AGROW>

                [J1 V1] = fast_track(I(:,idx),J0(:,idx),V0(:,idx));
                
                J(:,idx0) = J1(:,subidx);
                V(:,idx0) = V1(:,subidx);
            end
            V(:,end) = [];
            varargout = {J V};
        case 'fakelines'
            % Compute fake line images
            V = varargin{2};
            par = struct('vmax',ceil(max(abs(V(:)))));
            par = parameters(V,par);
            % (change sizes of J and V)
            %J = randn(par.nx3,par.nt); 
            J = repmat(sin(2*pi*(1:par.nx3)*.05)',1,par.nt);
            V = interp1(1:par.nx,V(:,1:end-1),linspace(1,par.nx,par.nx2));
            % (predict)
            for k=1:par.nt
                J(par.vmax2+(1:par.nx2),2:par.nt) = predict(J,V,par,false);
            end
            % (go back to coordinates as in I)
            J = J(par.vmax2+(1:par.nx2),:);
            varargout = {J};
        case 'par'
            % Default parameters
            par = defaultparameters;
            varargout{1} = par;
        otherwise
            error argument
    end
else
    % Estimation
    [I par] = deal(varargin{1:2});
    I = double(I);
    par = parameters(I,par);
    par.ha = [];
    % (initialization)
    if nargin<4
        [J0 V0 specialinit] = init(I,par);
    else
        [J0 V0] = deal(varargin{3:4});
        J0 = double(J0); V0 = double(V0);
        V0 = max(par.vmax(1),min(par.vmax(2),V0));
        specialinit = false;
    end
    % (init display)
    for i=3:nargin
        a = varargin{i};
        if isvector(a) && all(ishandle(a))
            par.ha = a;
        end
    end
    if isempty(par.ha)
        par.hf = 524;
        if ~ishandle(par.hf)
            figure(par.hf)
            set(par.hf,'numbertitle','off','name','track')
            par.ha(1) = axes('parent',524,'position',[.05 3/4+.04 .9 .20]);
            par.ha(2) = axes('parent',524,'position',[.05 2/4+.04 .9 .20]);
            par.ha(3) = axes('parent',524,'position',[.05 1/4+.04 .9 .20]);
            par.ha(4) = axes('parent',524,'position',[.05 .04 .43 .20]);
            par.ha(5) = axes('parent',524,'position',[.54 .04 .43 .20], ...
                'xlimmode','manual');
            par.hu = uicontrol('parent',524,'position',[10 10 120 15], ...
                'style','togglebutton','string','update display', ...
                'value',1);
        else
            par.ha = findobj(par.hf,'type','axes');
            par.ha = flipud(par.ha);
            cla(par.ha(5))
            par.hu = findobj(par.hf,'type','uicontrol');
            if get(par.hu,'value'), figure(par.hf), end
        end
    end
    par.him(1) = imagesc(I,'parent',par.ha(1)); grid(par.ha(1),'on');
    par.him(2) = imagesc(J0,'parent',par.ha(2)); grid(par.ha(2),'on');
    par.him(3) = imagesc(V0,'parent',par.ha(3)); grid(par.ha(3),'on');
    % (change sizes of J0 and V0)
    J0 = interp1(1:par.nx,J0, ...
        [linspace(1,1,par.vmax2) linspace(1,par.nx,par.nx2) linspace(par.nx,par.nx,par.vmax2)]);
    V0 = interp1(1:par.nx,V0(:,1:end-1),linspace(1,par.nx,par.nx2));
    % (estimation)
    [J V] = estimate(I,par,J0,V0,specialinit);
    % (go back to coordinates as in I)
    J = J(par.vmax2+(1:par.nx2),:);
    V = V(1:par.nx2,:); V = (V(:,[1 1:end])+V(:,[1:end end]))/2;
    % (final display)
    set(par.him(2),'cdata',J)
    set(par.him(3),'cdata',V)
    drawnow
    % (output)
    varargout = {J V};
end

%---------------
% MAIN FUNCTIONS
%---------------
function par = defaultparameters

% optimization parameters
par.kfit = 20;
par.kjt = 200;
par.kjx = 2;
par.kvt = 100;
par.kvx = 100;

% maximal speed (pixel/frame unit and subpixel/frame)
par.vmax = 5;

% initialization - if empty, estimated using fft
par.v0 = [];

%---
function par = parameters(I,par)
% input 'par' contains algorithms parameters
% output 'par' contains additional pre-computed values specific to the size
% of the image I

% min and max
if isscalar(par.vmax), par.vmax = par.vmax*[-1 1]; end

% number of additional spatial points ('frame')
par.vmax2 = ceil(max(abs(par.vmax))+1);

% general values
[par.nx par.nt] = size(I);
par.nx2 = 1 + (par.nx-1); % (sub-interpolation for V, K)
par.nx3 = par.nx2 + 2*par.vmax2; % (frames for J)

%---
function [J0 V0 specialinit] = init(I,par)

% J0: add 'frame'
J0 = double(I); %interp1(1:par.nx,I,[linspace(1,1,par.vmax2) linspace(1,par.nx,par.nx2) linspace(par.nx,par.nx,par.vmax2)]);

specialinit = false;
if ~isfield(par,'v0') || isempty(par.v0), par.v0 = 'fft'; end
if iscell(par.v0), disp('cell array for par.v0 is deprecated'), par.v0=par.v0{1}; end
if isnumeric(par.v0)
    v0 = double(par.v0);
else
    switch par.v0
        case 'fft'
            % initialization: estimate average speed using fft
            If = abs(fft2(I));
            If = If(2:end,2:end);
            Iff = filt2(If,2); m = max(Iff(:)); Iff(Iff<m/2)=0;
            fx = 1/par.nx * (mod((1:par.nx-1)+floor(par.nx/2),par.nx)'-floor(par.nx/2));
            ft = 1/par.nt * (mod((1:par.nt-1)+floor(par.nt/2),par.nt)-floor(par.nt/2));
            vv = repmat(ft,par.nx-1,1)./repmat(fx,1,par.nt-1);
            v0 = -mean(Iff(:).*vv(:)) / mean(Iff(:));
            disp(['Speed initialization at ' num2str(v0)])
        case 'est'
            specialinit = true;
            v0 = 0;
        otherwise
            error('unknown estimation flag ''%s''',par.v0)
    end
end
V0 = ones(par.nx,par.nt)*v0;

%---
function [J V] = estimate(I,par,J0,V0,specialinit)
% note that here, J0, V0, J and V do NOT have the same size as I

% normalization of I
fact = std(I(:));
I = I/fact;
J0 = J0/fact;

% more parameters
Jmax = inf(par.nx3,par.nt);
Vmin = par.vmax(1) * ones(par.nx2,par.nt-1);
xmin = interleave(-Jmax,Vmin,par);
Vmax = par.vmax(2) * ones(par.nx2,par.nt-1);
xmax = interleave(Jmax,Vmax,par);

% initialization
if specialinit
    % tries a number of initializations and choose the best
    vmin = par.vmax(1); vmax = par.vmax(2);
    [Es v0s nv] = testinit(vmin,vmax,.3,I,J0,V0,par,xmin,xmax);
    % repeat with a thinner grid and a narrower range
    [Es ord] = sort(Es);
    if ord(1)==1
        vmax = v0s(2);
    elseif ord(1)==nv
        vmin = v0s(nv-1);
    else %if all(abs(ord(2:3)-ord(1))==1)
        vmin = v0s(ord(1)-1);
        vmax = v0s(ord(1)+1);
    end
    [Es v0s] = testinit(vmin,vmax,.1,I,J0,V0,par,xmin,xmax);
    % choose the minimum without further check now
    [dum idx] = min(Es);
    v0 = v0s(idx);
    V0(:) = v0;
end
x0 = interleave(J0,V0,par);

% estimation
% (fmincon)
opt = optimoptions('fmincon','Algorithm','trust-region-reflective', ...
    'Gradobj','on','Hessian','on', ...
    'Display','none', ...
    'TolFun',2e-3, ...
    'MaxIter',100);
x = fmincon(@(x)energy(x,I,par),x0,[],[],[],[],xmin,xmax,[],opt);

% % (lsqnonlin)
% opt = optimset('LargeScale','on', ...
%     'Algorithm','levenberg-marquardt', ...
%     'Jacobian','on', ...
%     'Display','iter', ...
%     'MaxIter',300);
% x = lsqnonlin(@(x)myfun(x,I,par),x0,[],[],opt);

% terminate
[J V] = interleave(x,par);
J = J*fact;
if specialinit
    E = energy(x,I,par); v = mean(V(:));
    line(v,E,'parent',par.ha(5),'marker','*','color','k')
    drawnow
end

%---
function [Es speeds nv] = testinit(vmin,vmax,dvrel,I,J0,V0,par,xmin,xmax)
% tries several values for the initialization and returns the vector of
% resulting energies

% optimization parameters
opt = optimset('LargeScale','on', ...
    'Gradobj','on','Hessian','on', ...
    'Algorithm','trust-region-reflective', ...
    'Display','none', ...
    'TolFun',2e-3, ...
    'MaxIter',1);

% speeds to try
speeds = [];
dr = log(1+dvrel);
cut = 2;
if vmin<-cut
    rmin = log(-vmin);
    rmax = log(-min(vmax,-cut));
    r = linspace(rmin,rmax,1+ceil((rmin-rmax)/dr));
    speeds = [speeds -exp(r)];
end
if vmin<cut && vmax>-cut
    sidea = max(vmin,-cut); sideb = min(vmax,cut);
    speeds = [speeds linspace(sidea,sideb,1+ceil((sideb-sidea)/(dvrel*cut)))];
end
if vmax>cut
    rmin = log(max(vmin,cut));
    rmax = log(vmax);
    r = linspace(rmin,rmax,1+ceil((rmax-rmin)/dr));
    speeds = [speeds exp(r)];
end
speeds = unique(speeds);

% display: new plot (first grid) or plot on top of existing plot (second
% grid)?
if dvrel>=.3
    hl=[]; 
else
    hl=line(nan,nan,'parent',par.ha(5),'linewidth',2,'color','k'); 
end
% try all these elements
nv = length(speeds);
Es = zeros(1,nv);
for i=1:nv
    v0 = speeds(i);
    V0(:) = v0;
    x0 = interleave(J0,V0,par);
    x = fmincon(@(x)energy(x,I,par),x0,[],[],[],[],xmin,xmax,[],opt);
    Es(i) = energy(x,I,par);
    if isempty(hl)
        plot(speeds(1:i),Es(1:i),'parent',par.ha(5))
        set(par.ha(5),'xlim',speeds([1 nv]))
        if i>=4
            ylim = [min(Es(1:i)) max(Es(2:i-1))];
            ylim = ylim + [-.5 .5]*diff(ylim);
            set(par.ha(5),'ylim',ylim), 
        end
        drawnow
    else
        set(hl,'xdata',speeds(1:i),'ydata',Es(1:i)), drawnow
    end
end

%----------------------
% COMPUTATION FUNCTIONS
%----------------------
function [E G H] = energy(x,I0,par,dodetails)

if nargin<4, dodetails=false; end
if dodetails && nargout~=2, error('2 output arguments required when asking for details'), end

[J V] = interleave(x,par);
switch nargout
    case 1
        I = subsample(J,par);
        K = predict(J,V,par,false);
        L = predict(J,V,par,true);
        [Djt Djx Dvt Dvx] = differences(J,V,K,L,par);
    case 2
        [I dIdJ] = subsample(J,par);
        [K dKdJ dKdV] = predict(J,V,par,false);
        [L dLdJ dLdV] = predict(J,V,par,true);
        [Djt Djx Dvt Dvx dDjtdJ dDjxdJ dDvtdV dDvxdV] = differences(J,V,K,L,par,dKdJ,dLdJ);
        dDjtdV = [-dKdV; dLdV];
    case 3
        [I dIdJ] = subsample(J,par);
        [K dKdJ dKdV d2KdVdJ d2KdV2] = predict(J,V,par,false);  %#ok<ASGLU>
        [L dLdJ dLdV d2LdVdJ d2LdV2] = predict(J,V,par,true);  %#ok<ASGLU>
        [Djt Djx Dvt Dvx dDjtdJ dDjxdJ dDvtdV dDvxdV] = differences(J,V,K,L,par,dKdJ,dLdJ);
        dDjtdV = [-dKdV; dLdV];
end
Dfit = I-I0;

E = par.kfit * sum(Dfit(:).^2) ...  Dfit is the difference the estimated image of lines and the original image of lines
    + par.kjt * sum(Djt(:).^2) ...  Djt is the difference between successive columns (time instants) of the estimated image of lines, after translation (and re-interpolation) by the estimated speed
    + par.kjx * sum(Djx(:).^2) ...  Djx is the difference between successive rows (points) of the estimated image of lines
    + par.kvt * sum(Dvt(:).^2) ...  Dvt is the difference between successive columns (time instants) of the estimated speeds
    + par.kvx * sum(Dvx(:).^2);   % Dvx is the difference between successive rows (points) of the estimated speeds
if dodetails
    G = struct('fit',par.kfit*sum(Dfit(:).^2), ...
        'jt',par.kjt*sum(Djt(:).^2),'jx',par.kjx*sum(Djx(:).^2), ...
        'vt',par.kvt*sum(Dvt(:).^2),'vx',par.kvx*sum(Dvx(:).^2));
    return
end

if nargout>1
    dEdJ = 2*( par.kfit * (Dfit(:)'*dIdJ) ...
        + par.kjt * (Djt(:)'*dDjtdJ) ...
        + par.kjx * (Djx(:)'*dDjxdJ) );
    dEdV = 2*( par.kjt * (Djt(:)'*dDjtdV) ...
        + par.kvt * (Dvt(:)'*dDvtdV) ...
        + par.kvx * (Dvx(:)'*dDvxdV) );
    G = interleave(dEdJ,dEdV,par);
    if nargout>2
        % note that: 
        % e             = sum_k(f_k^2)
        % de/dx_i       = 2 * sum_k(f_k*df_k/dx_i)
        %               = 2 * f'*df/dx_i
        % (d2e/dx2)_ij  = 2 * sum_k(df_k/dx_i*df_k/dx_j + f_k*d2f_k/dx_idx_j)
        %               = 2 * (df/dx_i'*df/dx_j + f'*d2f/dx_idx_j)
        d2EdJ2 = 2*( par.kfit * (dIdJ'*dIdJ) ...
            + par.kjt * (dDjtdJ'*dDjtdJ) ...
            + par.kjx * (dDjxdJ'*dDjxdJ) );
        % we have d2K_k/dV_idJ_j non zero only if k=i
        % d2KdVdJ(i,j) = d2K_i/dV_idJ_j [see #1]
        % so d2EdVdJ(i,j) associated to K = Djt(i)*d2Djt_i/dV_idJ_j
        %                                 = - Djt(i)*d2KdVdJ(i,j)
        
        %tmp = d2KdVdJ;
        %[i j x] = find(tmp);
        %x = x .* Djt(i);
        %tmp = sparse(i,j,x,par.nx2*(par.nt-1),par.nx3*par.nt); 
        d2EdVdJ = 2 * ( par.kjt * dDjtdV'*dDjtdJ ...
            ); %- par.kjt * tmp ); MISSING BACKWARD DIRECTION!!
        
        % (now i stop explaining!!!)
        d2EdV2 = 2*( par.kjt * (dDjtdV'*dDjtdV) ...
            ... - par.kjt * spdiags(Djt(:).*d2KdV2,0,par.nx2*(par.nt-1),par.nx2*(par.nt-1)) MISSING BACKWARD DIRECTION!! ...
            + par.kvt * (dDvtdV'*dDvtdV) ...
            + par.kvx * (dDvxdV'*dDvxdV) );
        H = interleave(d2EdJ2,d2EdVdJ,d2EdV2,par);
        % display
        if get(par.hu,'value')
            switch length(par.ha)
                case 1
                    imagesc(I,'parent',par.ha)
                case 2
                    imagesc(I,'parent',par.ha(1))
                    imagesc(V(:,[1:end end]),'parent',par.ha(2))
                case 5
                    set(par.him(1),'cdata',I0)
                    set(par.him(2),'cdata',I)
                    set(par.him(3),'cdata',V(:,[1:end end]))
                    idx = round(par.nx2/3):round(par.nx2*2/3);
                    plot(mean(V(idx,:)),'parent',par.ha(4)), grid(par.ha(4),'on'), drawnow
            end
            drawnow
        end
    end
end

%---
function [I dIdJ] = subsample(J,par)
% note that this is a linear function

idx = 1 + par.vmax2 + (0:par.nx-1);
I = J(idx,:);

if nargout>1
    tmp = false(par.nx3,par.nt); 
    tmp(idx,:) = true;
    allidx = find(tmp);
    dIdJ = sparse(1:par.nx*par.nt,allidx,1,par.nx*par.nt,par.nx3*par.nt);
end

%---
function [K dKdJ dKdV d2KdVdJ d2KdV2] = predict(J,V,par,backward)
% backward is true [false] for backward [forward] prediction

if nargin<4, backward=true; end

% piece-wise polynomial interpolation: for x in [0 1],
% p(i+x) = [1 x x^2 x^3] * COEF * [f(i-1) f(i) f(i+1) f(i+2)]', with
% 1) p(i) = f(i)
% 2) p(i+1) = f(i+1),
% 3) p'(i) = (f(i+1)-f(i-1))/2
% 4) p'(i+1) = (f(i+2)-f(i))/2

COEF = [0 1 0 0; 0 -1 1 0; 0 0 0 0; 0 0 0 0]; % linear interpolations
% COEF = 1/2 * [0 2 0 0; -1 0 1 0; 2 -5 4 -1; -1 3 -3 1]; % cubic interpolations -> would necessit nx2 to be one more! 

if ~backward, V = -V; end
V0 = floor(V);
Vx = V-V0;

X{1} = ones(par.nx2,par.nt-1);
X{2} = Vx;
X{3} = Vx.*Vx;
X{4} = Vx.*X{3};

xidx = (1:par.nx2) + par.vmax2;
xidx = repmat(xidx',1,par.nt-1) + V0;
tidx = repmat((1:par.nt-1)+(backward),par.nx2,1);
allidx = xidx + par.nx3*(tidx-1);

%Ks{1} = J(allidx-1);
Ks{2} = J(allidx);
Ks{3} = J(allidx+1);
%Ks{4} = J(allidx+2);

K = zeros(par.nx2,par.nt-1);
for i=1:2 % 1:4
    for j=2:3 % 1:4
        a = COEF(i,j);
        if a, K = K + a*(X{i}.*Ks{j}); end
    end
end

% first order derivative?
if nargout>1
   % dKdJ(i,j) = d(K(i)) / d(J(j))
   dKdJ = sparse([],[],[],par.nx2*(par.nt-1),par.nx3*par.nt);
   for j=2:3 % 1:4
       tmp = zeros(par.nx2*(par.nt-1),1);
       for i=1:2 % 1:4
           a = COEF(i,j);
           if a, tmp = tmp + a*X{i}(:); end
       end
       dKdJ = dKdJ + sparse(1:par.nx2*(par.nt-1), allidx(:)+(j-2), tmp, ...
           par.nx2*(par.nt-1),par.nx3*par.nt);
   end
   
   % dKdV(i,i) = d(K(i)) / d(V(i))
   dKdV = zeros(par.nx2,par.nt-1);
   for i=2 % 2:4
       for j=2:3 % 1:4
           a = COEF(i,j);
           if a, dKdV = dKdV + ((i-1)*a)*(X{i-1}.*Ks{j}); end
       end
   end
   dKdV = spdiags(dKdV(:),0,par.nx2*(par.nt-1),par.nx2*(par.nt-1));
   if ~backward, dKdV = -dKdV; end
   
   % second order derivative?
   if nargout>2
       % d2KdVdJ(i,j) = d2(K(i)) / d(V(i)) d(J(j)) [#1]
       d2KdVdJ = sparse([],[],[],par.nx2*(par.nt-1),par.nx3*par.nt);
       for j=2:3 % 1:4
           tmp = zeros(par.nx2*(par.nt-1),1);
           for i=2 % 2:4
               a = COEF(i,j);
               if a, tmp = tmp + ((i-1)*a)*X{i-1}(:); end
           end
           d2KdVdJ = d2KdVdJ + sparse(1:par.nx2*(par.nt-1), allidx(:)+(j-2), tmp, ...
               par.nx2*(par.nt-1),par.nx3*par.nt);
       end
       if ~backward, d2KdVdJ = - d2KdVdJ; end
       
       % d2KdV2(i)  = d2(K(i)) / d(V(i))^2
       d2KdV2 = zeros(par.nx2,par.nt-1);
       for i=[] % 3:4
           for j=2:3 % 1:4
               a = COEF(i,j);
               if a, d2KdV2 = d2KdV2 + ((i-1)*(i-2)*a)*(X{i-2}.*Ks{j}); end
           end
       end
       d2KdV2 = d2KdV2(:);
   end
end

%---
function [Djt Djx Dvt Dvx dDjtdJ dDjxdJ dDvtdV dDvxdV] = differences(J,V,K,L,par,dKdJ,dLdJ)

Djt  = [J(par.vmax2+(1:par.nx2),2:par.nt)-K  L-J(par.vmax2+(1:par.nx2),1:par.nt-1)];
Djx  = J(2:end,:)-J(1:end-1,:);
Dvt  = V(:,2:end)-V(:,1:end-1);
Dvx  = V(2:end,:)-V(1:end-1,:);

if nargout>4
    allidx1 = repmat(par.vmax2+(1:par.nx2)',1,par.nt-1) ...
        + par.nx3*repmat(0:par.nt-2,par.nx2,1);
    allidx2 = repmat((1:par.nx2)',1,par.nt-1) ...
        + par.nx2*repmat(0:par.nt-2,par.nx2,1);
    dDjtdJ1 = sparse(allidx2,allidx1+par.nx3,1,par.nx2*(par.nt-1),par.nx3*par.nt) ...
        - dKdJ;
    dDjtdJ2 = dLdJ ...
        - sparse(allidx2,allidx1,1,par.nx2*(par.nt-1),par.nx3*par.nt);
    dDjtdJ = [dDjtdJ1; dDjtdJ2];

    allidx1 = repmat((1:par.nx3-1)',1,par.nt) ...
        + par.nx3*repmat(0:par.nt-1,par.nx3-1,1);
    allidx2 = repmat((1:par.nx3-1)',1,par.nt) ...
        + (par.nx3-1)*repmat(0:par.nt-1,par.nx3-1,1);
    dDjxdJ = sparse(allidx2,allidx1+1,1,(par.nx3-1)*par.nt,par.nx3*par.nt) ...
        - sparse(allidx2,allidx1,1,(par.nx3-1)*par.nt,par.nx3*par.nt);

    allidx1 = repmat((1:par.nx2)',1,par.nt-2) ...
        + par.nx2*repmat(0:par.nt-3,par.nx2,1);
    allidx2 = repmat((1:par.nx2)',1,par.nt-2) ...
        + par.nx2*repmat(0:par.nt-3,par.nx2,1);
    dDvtdV = sparse(allidx2,allidx1+par.nx2,1,par.nx2*(par.nt-2),par.nx2*(par.nt-1)) ...
        - sparse(allidx2,allidx1,1,par.nx2*(par.nt-2),par.nx2*(par.nt-1));

    allidx1 = repmat((1:par.nx2-1)',1,par.nt-1) ...
        + par.nx2*repmat(0:par.nt-2,par.nx2-1,1);
    allidx2 = repmat((1:par.nx2-1)',1,par.nt-1) ...
        + (par.nx2-1)*repmat(0:par.nt-2,par.nx2-1,1);
    dDvxdV = sparse(allidx2,allidx1+1,1,(par.nx2-1)*(par.nt-1),par.nx2*(par.nt-1)) ...
        - sparse(allidx2,allidx1,1,(par.nx2-1)*(par.nt-1),par.nx2*(par.nt-1));
end

%---
function varargout = interleave(varargin)
% function x = interleave(J,V,par);
% function [J V] = interleave(x,par);
% function G = interleave(dEdJ,dEdV,par);
% function H = interleave(d2EdJ2,d2EdVdJ,d2EdV2,par);
% function Jac = interleave(dfdJ,dfdV,par);

par = varargin{end};

nt = par.nt;
nx2 = par.nx2;
nx3 = par.nx3;

% technical notes:
% - in order for the algorithm to initially set a large trust region,
%   we 'reduce' the norm of x by a factor (for example, 10)

FACT = 10;

if nargout==2
    % function [J V] = interleave(x,par);
    x = varargin{1}*FACT;
    J = reshape(x(1:nx3*nt),nx3,nt);
    V = reshape(x(nx3*nt+(1:nx2*(nt-1))),nx2,nt-1);
    varargout = {J V};
elseif nargin==4
    % function H = interleave(d2EdJ2,d2EdVdJ,d2EdV2,par);
    [d2EdJ2 d2EdVdJ d2EdV2] = deal(varargin{1:3});
    H = [d2EdJ2 d2EdVdJ'; d2EdVdJ d2EdV2];
    varargout = {H*(FACT*FACT)};
elseif isvector(varargin{1})
    % function G = interleave(dEdJ,dEdV,par);
    [dEdJ dEdV] = deal(varargin{1:2});
    G = [dEdJ dEdV];
    varargout = {G*FACT};
elseif size(varargin{1},2)==nx3*nt
    % function J = interleave(dfdJ,dfdV,par);
    [dfdJ dfdV] = deal(varargin{1:2});
    Jac = [dfdJ dfdV];
    varargout = {Jac*FACT};
else
    % function x = interleave(J,V,par);
    [J V] = deal(varargin{1:2});
    x = [J(:); V(:)];
    varargout = {x/FACT};
end


