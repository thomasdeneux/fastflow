function [x xm] = fast_radon(I,par)
% function par = fast_radon('par') [get default parameters]
% function [V v] = fast_radon(I,par)
%---
% v is an average value computed by using all the information in I; it
% should be a better estimation than mean(V(:))

if nargin==0, help fast_radon, return, end

% default parameters
if nargin==1 && strcmp(I,'par')
    x = defaultparameters;
    return
end

% estimation
if nargin<2, par = defaultparameters; end
if nargout<2
    x = estimate(I,par);
% elseif isscalar(par.delta)
%     x = estimate_old(I,par);
else
    x = estimate_old(I,par);
end

%---
function par = defaultparameters

% % resample
% par.rsampx = 1;
% par.rsampt = 1;

% box size
par.delta = [50 50];

% box shifts
par.shiftx = 1/2;
par.shiftt = 1/4;

% angles to try
par.dvrel = .01; % in degrees
par.vmax = 15;

% % display 
% par.display = true;


%---
function V = estimate(I,par)
% CONCEPT
% The Radon transform can be computed on any rectangle, not only
% a square; this is achieved by sliding a square inside the rectangle and
% averaging the successive Radon transforms for all square positions!
% Or rather, the variance of the Radon transforms.
% In the comments, this square is called "the square" and its size is
% delta, whereas the rectangle is called "the Radon box" and its size is
% [deltax deltat].
% The parameter par.delta specified by user stands for [deltax deltat] if
% it contains 2 values, or deltat if it is scalar (then deltax = nx). The
% parameters par.shiftx and par.shiftt specify the minimal shfits, i.e. the
% required resolution.

[nx nt] = size(I);

% angles
if isscalar(par.vmax), par.vmax=par.vmax*[-1 1]; end
if isfield(par,'dvrel')
    % values of speed: increase by dvrel/10 for abs(v)<1/10, and by v*dvrel
    % for abs(v)>1/10
    doang = false;
    vmin = par.vmax(1); vmax = par.vmax(2);
    speeds = [];
    dr = log(1+par.dvrel);
    cut = .7;
    if vmin<-cut
        rmin = log(-vmin);
        rmax = log(-min(vmax,-cut));
        r = linspace(rmin,rmax,1+ceil((rmin-rmax)/dr));
        speeds = [speeds -exp(r)];
    end
    if vmin<cut && vmax>-cut
        speeds = [speeds max(vmin,-cut):par.dvrel*cut:min(vmax,cut)];
    end
    if vmax>cut
        rmin = log(max(vmin,cut));
        rmax = log(vmax);
        r = linspace(rmin,rmax,1+ceil((rmax-rmin)/dr));
        speeds = [speeds exp(r)];
    end
    speeds = unique(speeds);
    angles = mod(acotd(speeds),180);
elseif isfield(par,'dtheta')
    doang = true;
    angmax = fn_round(mod(acotd(par.vmax),180),par.dtheta);
    angles = angmax(2):par.dtheta:angmax(1); % ok since mod(acotd(x),180) is decreasing
    speeds = cotd(angles);
else 
    error('parameter missing')
end
nv = length(speeds);

% box parameters specified by user
if isscalar(par.delta), par.delta = [nx par.delta]; end 
deltax = min(par.delta(1),nx);
deltat = min(par.delta(2),nt);
shiftx = par.shiftx; if shiftx<1, shiftx = floor(deltax*shiftx); end
shiftt = par.shiftt; if shiftt<1, shiftt = floor(deltat*shiftt); end

% size of the sliding square
delta  = min(deltax,deltat);
% (it is ridiculous to have a very small shift, therefore, increase the
% square side to nx or nt if it is already close to it)
if delta > min(nx,nt)/1.5
    delta  = min(nx,nt);
    % the following keeps 'delta = min(deltax,deltat)' true
    deltax = max(delta,deltax); 
    deltat = max(delta,deltat);
end

% how to shift the square
% (shift of the square required to cover the Radon box)
if delta<deltax
    shiftx = min(shiftx,delta/2);
elseif delta<deltat
    shiftt = min(shiftt,delta/2);
end
% (adjust to have an integer number of squares spanning the whole image)
nsquarex  = 1+ceil((nx-delta)/shiftx); 
nsquaret  = 1+ceil((nt-delta)/shiftt); 
shiftx = (nx-delta)/(nsquarex-1);
shiftt = (nt-delta)/(nsquaret-1);
startx = round(linspace(0,nx-delta,nsquarex));
startt = round(linspace(0,nt-delta,nsquaret));

% how many squares inside a Radon box, and how many Radon boxes
nsqx = 1+round((deltax-delta)/shiftx); if nsquarex==1, nsqx=1; end
nsqt = 1+round((deltat-delta)/shiftt); if nsquaret==1, nsqt=1; end
nboxx = nsquarex-(nsqx-1);
nboxt = nsquaret-(nsqt-1);
    
% other parameters
mask = fspecial('gaussian',delta,delta/4);
[dum Xp] = radon(zeros(delta),angles); rkeep = (abs(Xp)<delta/2);

% prepare display
hf = 710;
if ~ishandle(hf)
    figure(hf)
    set(hf,'numbertitle','off','name','radon')
    hudisplay = uicontrol('style','togglebutton','string','update display', ...
        'value',1, ...
        'position',[10 10 120 20]);
else
    hudisplay = findobj(hf,'type','uicontrol');
    if ~isscalar(hudisplay), error('close window and start again'), end
    if get(hudisplay,'value'), figure(hf), end
end
ha = subplot(2,1,1,'parent',hf);
hb = subplot(2,2,3,'parent',hf); 
hc = subplot(2,2,4,'parent',hf); 
drawnow
hrec = [];

% Slide the square!
RVS = zeros(nv,nsquarex,nsquaret);
for ix=1:nsquarex
    idxx = startx(ix)+(1:delta);
    for it=1:nsquaret
        % cut box
        idxt = startt(it)+(1:delta);
        box = I(idxx,idxt);
        % estimate
        box = box.*mask;
        box = box - (sum(box(:))/sum(mask(:))).*mask;
        R = radon(box,angles);
        curvar = var(R(rkeep,:)); 
        RVS(:,ix,it) = curvar;
        [dum imax] = max(curvar);
        % display
        if ishandle(hf) && get(hudisplay,'value')
            if isempty(hrec)
                [hrec hib hlb hic hlc] = firstdisplay(I,delta,angles,doang,speeds,ha,hb,hc);
            end
            set(hrec,'position',[startt(it)+.5 startx(ix)+.5 delta delta])
            set(hib,'cdata',box)
            theta = angles(imax);
            set(hlb,'xdata',delta/2*(1+[-1 1]*sind(theta)), ...
                'ydata',delta/2*(1+[-1 1]*cosd(theta)))
            set(hic,'cdata',R(rkeep,:))
            if doang
                set(hlc,'xdata',[1 1]*theta)
            else
                set(hlc,'xdata',[1 1]*imax)
            end
            drawnow
        end
    end
end
delete(hrec)

% Now, slide the Radon box!
RVB = RVS(:,1:nboxx,1:nboxt);
for i=2:nsqx, RVB = RVB + RVS(:,i:nboxx+i-1,:); end
for i=2:nsqt, RVB = RVB + RVS(:,:,i:nboxt+i-1); end

% Finally collect the velocity estimates
V = zeros(nx,nt);
startt = round(linspace(0,nt-deltat,nboxt));
for ix=1:nboxx
    idxx = floor(nx*(ix-1)/nboxx)+1:floor(nx*ix/nboxx);
    if ix==1
        idxx = 1:deltax;
    else
        idxx = startx(ix)+(1+round(deltax/2-shiftx/2):deltax);
    end
    for it=1:nboxt
        % fill result
        if it==1
            idxt = 1:deltat;
        else
            idxt = startt(it)+(1+round(deltat/2-shiftt/2):deltat); 
        end
        [dum imax] = max(RVB(:,ix,it));
        V(idxx,idxt) = speeds(imax);
    end 
end

%---
function [hrec hib hlb hic hlc] = firstdisplay(I,delta,angles,doang,speeds,ha,hb,hc)

imagesc(I,'parent',ha);
hrec = rectangle('parent',ha,'edgecolor','b','erasemode','xor');
box = zeros(delta);
hib=imagesc(box,'parent',hb);
hlb = line(0,0,'parent',hb,'erasemode','xor');
[R Xp] = radon(box,angles);
rkeep = (abs(Xp)<delta/2);
if doang
    hic=imagesc(angles,Xp(rkeep),R(rkeep,:),'parent',hc);
else
    hic=imagesc(1:length(speeds),Xp(rkeep),R(rkeep,:),'parent',hc);
    ticks = [-20 -15 -10 -7 -5:-1 -.5 -.1 .1 .5 1:5 7 10 15 20];
    ticks = ticks(ticks>speeds(1) & ticks<speeds(end));
    idx = zeros(1,length(ticks));
    for i=1:length(ticks), idx(i) = find(speeds>=ticks(i),1); end
    set(hc,'xtick',idx,'xticklabel',cellstr(num2str(ticks')))
end
hlc = line([0 0],Xp([1 end]),'parent',hc,'erasemode','xor');
axis(hc,'normal')


%---
function [V v] = estimate_old(I,par)

% if par.rsampx~=1 || par.rsampt~=1
%     error('resampling not implemented yet')
% end
dov = (nargout>=2);

% boxes and where to put them
[nx nt] = size(I);

% box parameters
% (delta smaller than image size, but also increase if close to x size)
delta = min(min(nx,nt),par.delta);
if nx<1.5*delta, delta = min(nx,nt); end % increase a bit delta to make nsquarex=1
% (how to shift the box)
shiftx = par.shiftx; shiftt = par.shiftt;
if shiftx<1, shiftx = delta*shiftx; end
if shiftt<1, shiftt = delta*shiftt; end
nsquarex  = 1+ceil((nx-delta)/shiftx);
nsquaret  = 1+ceil((nt-delta)/shiftt);
%shiftx = (nx-delta)/(nsquarex-1);
shiftt = (nt-delta)/(nsquaret-1);
startx = round(linspace(0,nx-delta,nsquarex));
startt = round(linspace(0,nt-delta,nsquaret));

% angles
if isscalar(par.vmax), par.vmax=par.vmax*[-1 1]; end
if isfield(par,'dvrel')
    % values of speed: increase by dvrel/10 for abs(v)<1/10, and by v*dvrel
    % for abs(v)>1/10
    doang = false;
    vmin = par.vmax(1); vmax = par.vmax(2);
    speeds = [];
    dr = log(1+par.dvrel);
    cut = .7;
    if vmin<-cut
        rmin = log(-vmin);
        rmax = log(-min(vmax,-cut));
        r = linspace(rmin,rmax,1+ceil((rmin-rmax)/dr));
        speeds = [speeds -exp(r)];
    end
    if vmin<cut && vmax>-cut
        speeds = [speeds max(vmin,-cut):par.dvrel*cut:min(vmax,cut)];
    end
    if vmax>cut
        rmin = log(max(vmin,cut));
        rmax = log(vmax);
        r = linspace(rmin,rmax,1+ceil((rmax-rmin)/dr));
        speeds = [speeds exp(r)];
    end
    speeds = unique(speeds);
    angles = mod(acotd(speeds),180);
elseif isfield(par,'dtheta')
    doang = true;
    angmax = fn_round(mod(acotd(par.vmax),180),par.dtheta);
    angles = angmax(2):par.dtheta:angmax(1); % ok since mod(acotd(x),180) is decreasing
    speeds = cotd(angles);
else 
    error('parameter missing')
end
% other parameters
mask = fspecial('gaussian',delta,delta/4);
% estimatio of average speed
if dov, mvar = 0; end

% display
hf = 710;
if ~ishandle(hf)
    figure(hf)
    set(hf,'numbertitle','off','name','radon')
    hudisplay = uicontrol('style','togglebutton','string','update display', ...
        'value',1, ...
        'position',[10 10 120 20]);
else
    hudisplay = findobj(hf,'type','uicontrol');
    if ~isscalar(hudisplay), error('close window and start again'), end
    if get(hudisplay,'value'), figure(hf), end
end
ha = subplot(2,1,1,'parent',hf);
imagesc(I,'parent',ha);
hrec = rectangle('parent',ha,'edgecolor','b','erasemode','xor');
hb=subplot(2,2,3,'parent',hf); 
box = zeros(delta);
hib=imagesc(box,'parent',hb);
hlb = line(0,0,'parent',hb,'erasemode','xor');
hc=subplot(2,2,4,'parent',hf); 
[R Xp] = radon(box,angles);
rkeep = (abs(Xp)<delta/2);
if doang
    hic=imagesc(angles,Xp(rkeep),R(rkeep,:),'parent',hc);
else
    hic=imagesc(1:length(speeds),Xp(rkeep),R(rkeep,:),'parent',hc);
    ticks = [-20 -15 -10 -7 -5:-1 -.5 -.1 .1 .5 1:5 7 10 15 20];
    ticks = ticks(ticks>vmin & ticks<vmax);
    idx = zeros(1,length(ticks));
    for i=1:length(ticks), idx(i) = find(speeds>=ticks(i),1); end
    set(hc,'xtick',idx,'xticklabel',cellstr(num2str(ticks')))
end
hlc = line([0 0],Xp([1 end]),'parent',hc,'erasemode','xor');
axis(hc,'normal')

% go!
V = zeros(nx,nt);
for ix=1:nsquarex
    idxx = startx(ix)+(1:delta);
    %if ix==1, idxx_in=idxx; else idxx_in=startx(ix)+(1+round(delta/2-shiftx/2):delta); end
    idxx_in = floor(nx*(ix-1)/nsquarex)+1:floor(nx*ix/nsquarex);
    for it=1:nsquaret
        % cut box
        idxt = startt(it)+(1:delta);
        box = I(idxx,idxt);
        % estimate
        box = box.*mask;
        box = box - (sum(box(:))/sum(mask(:))).*mask;
        R = radon(box,angles);
        curvar = var(R(rkeep,:));
        if dov, mvar = mvar+curvar; end
        [dum imax] = max(curvar);
        v = speeds(imax);
        % fill result
        if it==1, idxt_in=idxt; else idxt_in=startt(it)+(1+round(delta/2-shiftt/2):delta); end
        V(idxx_in,idxt_in) = v;
        % display
        if ishandle(hf) && get(hudisplay,'value')
            set(hrec,'position',[startt(it)+.5 startx(ix)+.5 delta delta])
            set(hib,'cdata',box)
            theta = angles(imax);
            set(hlb,'xdata',delta/2*(1+[-1 1]*sind(theta)), ...
                'ydata',delta/2*(1+[-1 1]*cosd(theta)))
            set(hic,'cdata',R(rkeep,:))
            if doang
                set(hlc,'xdata',[1 1]*theta)
            else
                set(hlc,'xdata',[1 1]*imax)
            end
        end
        drawnow
    end
end
delete(hrec)

% estimation of average speed
if dov, [dum imax] = max(mvar); v = speeds(imax); end


