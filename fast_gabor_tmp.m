function x = fast_gabor(varargin)
% [anglesc,angl,coef]=fast_gabor(data,par)
% par = fast_gabor('par')
%
%INPUTS:
%data : data to analyse
%aff: 0 for no graphical results. 1: for drawing a graph. 2: for drawing at each step
%filter_struct: Gabor filters
%elim90=1 : flag pour eliminer la moiti� des orientations en function de la moyenne sur tout les X et les t.
%elim90=2 :elimine la moiti� des orientations avec un test en chaque point.
%elim90=0 :ne fait rien.
%
%OUTPUTS:
%angl(t,x): results given as the angle in function of the time (media from 4 first scales).
%anglesc(t,x,sc): results scale by scale
%
%--- comment by Thomas Deneux ---
%   coef(:,:,sc,or)     dot product result with Gabor coefficient with
%                       scale sc and orientation or
%   anglesc_syl{sc}     average orientation index (using coef for
%                       weighting) at scale sc
%                       [note: use 'orv' to convert from orientation index
%                       to angle]
%   speedsc_syl{sc}     average tangente of the angle at scale sc
%   poids_syl{sc}       total weight in each pixel (sum of coef)
%
%   angletg_syl{sc}     average orientation index - correction by the 'bin
%                       width' (each bin covers the same DeltaAngle; as a
%                       result bins closer to 90deg cover a larger DeltaTg
%                       than bins closer to 0deg; an additional bin width
%                       correction factor is used then) - at scale sc
%                       [note that this bin correction is awkward when
%                       averaging angles, but meaningful when averaging
%                       tangentes, as just below]
%   speedtg_syl{sc}     average tangente of the angle - using bin width
%                       correction - at scale sc
%   poidstg_syl{sc}     total weight in each pixel - using bin width
%                       correction
%
%   anglesc_orig        average orientation over different scales
%   speedsc             average tangent of angles over different scales
%   poids_orig          total weight in each pixel over orientations and
%                       scales
%
%   angletg             average orientation using bin width correction over
%                       different scales
%   speedtg             average tangent of angles using bin width
%                       correction over different scales
%   poidstg             total weight in each pixel using bin width
%                       correction over orientations and scales
%
%   anglesc             average orientation over different scales -
%                       additional algorithm used to "fill the holes" were
%                       no estimation is reliable
%   poids               total weight in each pixel over orientations and
%                       scales, after this filling algorithm
%
%   bestang             orientation corresponding to the maximum coef (over
%                       orientations and scales)
%--- end of comment ---
%
% Sylvain Fischer, avril 2006, revised July 2006
% Sylvain Takerkart, march 2007
% Thomas Deneux, 2008, 2009

if nargin==1 && isequal(varargin{1},'par')
    % Default parameters
    x = defaultparameters;
else
    % Estimate (filters reconstructed each time because of memory
    % concerns...)
    [data par] = deal(varargin{:});
    [ni nj] = size(data);
    angles = par.min_angle + ...
        (0:par.nbr_orient-1)/(par.nbr_orient-1)*(par.max_angle-par.min_angle);
    nper = length(par.periods); %#ok<NASGU>
    next = length(par.widthfact);
    nlen = length(par.length);
    periods = kron(par.periods(:)',ones(1,next*nlen));
    width = fn_mult(par.periods(:)',par.widthfact(:));
    extents = kron(width(:)',ones(1,nlen));
    ellips = fn_mult(width(:)',1./par.length(:));
    ellips = ellips(:)';
    filter_size = pow2(floor(log2(max(ni,nj))) + 1);
    aff = par.gabordisplay;
    x = gabor_estimate(data,filter_size,periods,extents,ellips,angles,aff);
end

%---
function par = defaultparameters

par.nbr_orient = 89;
par.min_angle = 1;
par.max_angle = 89;
par.periods = 4 * power(1.3,[0 1 2]);
par.widthfact = [.5 1 2];
par.length = 8;
par.useresult = 'anglesc';
par.gabordisplay = 0;

%---
function res = gabor_estimate(data,filter_size,periods,extents,ellips,angles,aff)

%-- GABOR BASE PARAMETERS --% 

% Gabor parameters
nbr_period = length(periods);
nbr_extent = length(extents);
nbr_ellips = length(ellips);
nbr_scales = max([nbr_period nbr_extent nbr_ellips]);
if ~all([nbr_period nbr_extent nbr_ellips]==nbr_scales)
    error('number of periods, extents and ellipcities should all be the same or be 1')
end
nbr_orient=length(angles);
% M=(nbr_scales*nbr_orient)*filter_size^2; % overcompleteness
min_angle = angles(1);
max_angle = angles(end);

% coordinates
% [fY fX] = meshgrid(1:filter_size,1:filter_size);
[Y X] = meshgrid(-(filter_size)/2:(filter_size)/2-1,-(filter_size)/2:(filter_size)/2-1);

%-- ESTIMATION PARAMETERS --%

seuilPoids=10 ; %weight threshold to consider the weigthed sum as reliable

img1200b=filtx2(data,15,'hm');
tmax=size(img1200b,1) ;
xmax=size(img1200b,2) ;

wdwatt=8 ; %attenuation en x o t: taille de la zone de transition (i.e. du raised cosine)
% wdw=filter_size;    %attenuation 8 pixels bords donnees
att=(.5-.5*cos([1:wdwatt]'/wdwatt*pi)); %#ok<*NBRAK>
img1200b(1:wdwatt,:)=img1200b(1:wdwatt,:).*(att*ones(1,xmax)) ;
img1200b(tmax+1-wdwatt:tmax,:)=img1200b(tmax+1-wdwatt:tmax,:).*(att(wdwatt:-1:1)*ones(1,xmax)) ;
szdata=size(img1200b) ;
att=(.5-.5*cos([1:wdwatt]/wdwatt*pi)).^1 ;
img1200b(:,1:wdwatt)=img1200b(:,1:wdwatt).*(ones(tmax,1)*att) ;
img1200b(:,szdata(2)+1-wdwatt:szdata(2))=img1200b(:,szdata(2)+1-wdwatt:szdata(2)).*(ones(tmax,1)*att(wdwatt:-1:1)) ;
myim = zeros(filter_size);
tstart = floor((filter_size - tmax)/2);
tend = tstart + tmax - 1;
xstart = floor((filter_size - xmax)/2);
xend = xstart + xmax - 1;
myim(tstart:tend, xstart:xend) = img1200b;
FI=fft2(myim);
Ic=zeros(tmax,xmax,nbr_scales,nbr_orient);

%-- ESTIMATION LOOP --%

% determine in which quarter we are
If = abs(FI);
take = filter_size/2;
balance = mean(mean(If(2:take,2:take))) - mean(mean(If(2:take,take+2:2*take)));
if balance>0
    % "other" quarter
    disp('detected negative speed')
    angles = angles + 90;
else
    disp('detected positive speed')
end

fn_progress('Gabor weight',nbr_scales*nbr_orient)
i_total = 0;
for i_orient=1:nbr_orient
    theta = angles(i_orient);
    U = X*sind(theta) - Y*cosd(theta);
    V = X*cosd(theta) + Y*sind(theta);
    for i_scale=1:nbr_scales
        % Compute filter
        i_total = i_total+1;
        fn_progress(i_total)
        lambda = periods(i_scale);
        sigma = extents(i_scale);
        gamma = ellips(i_scale);
        ATTENUATION = exp(-(U.^2 + (V*gamma).^2) / sigma^2);
        PHASE = exp(2*pi*1i*U/lambda);
        GABOR = ATTENUATION .* PHASE;
        N_F=sum(sum(GABOR.*conj(GABOR)))/filter_size/filter_size;
        Gabor = fftshift(GABOR)/sqrt(N_F);
        Fk = fft2(Gabor);

        % Estimate
        N_win=512;
        tmpim=ifft2(FI.*conj(Fk))/N_win;
        Ic(:,:,i_scale,i_orient)=tmpim(tstart:tend, xstart:xend);
    end;
end;

coef=abs(Ic).^2 ;
tam=size(coef) ;
tmax=size(coef,1) ;
xmax=size(coef,2) ;

coef0=coef; % save the energies before ponderation
meannrj=mean(mean(mean(coef(1+wdwatt:tmax-wdwatt,1+wdwatt:xmax-wdwatt,:,:),1),2),4) ;
for sc=1:tam(3), %ponderar por la energia media en el canal
    coef(:,:,sc,:)=coef(:,:,sc,:)/meannrj(sc) ;
end
% meanpeaknrj=mean(mean(max(coef(1+wdwatt:tmax-wdwatt,1+wdwatt:xmax-wdwatt,:,:),[],4),1),2) ;
%meanpeaknrj=meanpeaknrj(1:tam(3)-1) ;
%mediannrj=mean(coef,4) ; %median for each position/scale across orientations

coef=coef-1 ; %baseline suppression
coef(coef<0)=0; %baseline suppression

orv = angles;

% compute bin width (as suggested per ivo)
bw(1) = (tand((orv(2)+orv(1))/2) - tand(orv(1)) ) * 2;
for or=2:tam(4)-1
    bw(or) = tand((orv(or+1)+orv(or))/2) - tand((orv(or)+orv(or-1))/2);
end;
bw(tam(4)) = ( tand(orv(tam(4))) - tand( (orv(tam(4))+orv(tam(4)-1)) / 2 ) ) * 2;


poids = zeros(tam(1:2));
anglesc = zeros(tam(1:2));
% speedsc = zeros(tam(1:2));
% poidstg = zeros(tam(1:2));
% angletg = zeros(tam(1:2));
% speedtg = zeros(tam(1:2));
for sc=1:tam(3)
    anglesc_syl{sc} = zeros(tam(1:2)); %#ok<*AGROW>
    %     speedsc_syl{sc} = zeros(tam(1:2));
    poids_syl{sc} =  zeros(tam(1:2));
    %     angletg_syl{sc} = zeros(tam(1:2));
    %     speedtg_syl{sc} = zeros(tam(1:2));
    %     poidstg_syl{sc} =  zeros(tam(1:2));
    for or=1:tam(4)
        anglesc_syl{sc} = anglesc_syl{sc} + coef(:,:,sc,or)*orv(or);
        %         speedsc_syl{sc} = speedsc_syl{sc} + coef(:,:,sc,or)*tand(orv(or));
        poids_syl{sc} = poids_syl{sc} + coef(:,:,sc,or);
        %         angletg_syl{sc} = angletg_syl{sc} + coef(:,:,sc,or)*orv(or)*bw(or);
        %         speedtg_syl{sc} = speedtg_syl{sc} + coef(:,:,sc,or)*tand(orv(or))*bw(or);
        %         poidstg_syl{sc} = poidstg_syl{sc} + coef(:,:,sc,or)*bw(or);
    end
    anglesc = anglesc + anglesc_syl{sc};
    %     speedsc = speedsc + speedsc_syl{sc};
    anglesc_syl{sc} = anglesc_syl{sc}./(poids_syl{sc}+1e-33);
    %     speedsc_syl{sc} = speedsc_syl{sc}./(poids_syl{sc}+1e-33);
    poids = poids + poids_syl{sc};
    
    %     angletg = angletg + angletg_syl{sc};
    %     speedtg = speedtg + speedtg_syl{sc};
    %     angletg_syl{sc} = angletg_syl{sc}./(poidstg_syl{sc}+1e-33);
    %     speedtg_syl{sc} = speedtg_syl{sc}./(poidstg_syl{sc}+1e-33);
    %     poidstg = poidstg + poidstg_syl{sc};
    
end
disp(['Orientation not defined in ',int2str(sum(sum(poids<=0))), ' points on ',int2str(prod(tam(1:2)))]) ; drawnow
disp(['Orientation unreliably defined in ',int2str(sum(sum(poids<seuilPoids))), ' points on ',int2str(prod(tam(1:2)))]) ; drawnow
anglesc=anglesc./(poids+1e-33) ; %weighted sum 
% angletg=angletg./(poidstg+1e-33) ; %weighted sum 
% speedsc=speedsc./(poids+1e-33) ; %weighted sum
% speedtg=speedtg./(poidstg+1e-33) ; %weighted sum 

poids_orig = poids;
anglesc_orig = anglesc;


%%%%%%%%%%%%%%  averaging for the less defined coef.
disp('Propagating neighboring values to pixels for which the angle is undefined or unreliable:'), drawnow
% cvK= [1    10    45   120   210   252   210   120    45    10     1]/1024 ; 
cvK= [1     8    28    56    70    56    28     8     1]/256 ; %[1 6 15 20 15 6 1]/64 ; %[1 4 6 4 1]
it=0;
while it<30 || ~all(poids(:))
    it = it+1;
    anglescM=conv2(anglesc.*poids,cvK'*cvK,'same') ; %convolucion con los vecinos
    poidsM=conv2(poids,cvK'*cvK,'same') ; %convolucion con los vecinos
    anglescM=anglescM./(poidsM+1e-33) ; %recalculo por pesos
    anglesc(poids<seuilPoids)=anglescM(poids<seuilPoids) ;%reemplzar los mal definidos
    poids(poids<seuilPoids)=poidsM(poids<seuilPoids) ; %reemplzar los mal definidos
    maldef=sum(sum(poids<seuilPoids)) ; %number of points unreliably defined
    if maldef==0, break ; %every point is already  accuretely definied
    end
end
disp([int2str(it),' : orientation unreliably defined in ',int2str(maldef), ' points on ',int2str(prod(tam(1:2)))]) ;
    

res.poids = poids;
res.poids_orig = poids_orig;
res.anglesc = anglesc;
res.anglesc_orig = anglesc_orig;
res.poids_syl = poids_syl;
res.anglesc_syl = anglesc_syl;
% res.speedsc = speedsc;
% res.speedsc_syl = speedsc_syl;

% res.poidstg = poidstg;
% res.angletg = angletg;
% res.poidstg_syl = poidstg_syl;
% res.angletg_syl = angletg_syl;
% res.speedtg = speedtg;
% res.speedtg_syl = speedtg_syl;

res.bw = bw; 

% % save Ic and coef
% disp('saving Ic and coef')
% res.Ic = Ic;
% res.coef = coef;

% find for each pixel the best filter
disp('saving best angle'), drawnow
[nt nx nsc ndir] = size(coef0);
% n = nt*nx;
c = coef0;
c = reshape(c,[nt*nx nsc*ndir]);
[dum idx] = max(c,[],2);
idx = reshape(idx,nt,nx);
% isc = 1+mod((idx-1),nsc);
idir = 1+floor((idx-1)/nsc);
angmin = min_angle;
angmax = max_angle;
ang = angmin + (0:(nbr_orient-1))/(nbr_orient-1) * (angmax-angmin);
% ind = 1:nt*nx;
% f = find(idir(:)~=1);
% recons = zeros(nt,nx);
% recons(ind(f)) = Ic(ind(f) + n*(idx(ind(f))-1));
% res.idx = idx;
% res.isc = isc;
% res.idir = idir;
res.bestang = ang(idir);
% res.recons = recons;

%---
function FGabor = test_base(par,ni,nj)

angles = par.min_angle + ...
    (0:par.nbr_orient-1)/(par.nbr_orient-1)*(par.max_angle-par.min_angle);
nper = length(par.periods); 
next = length(par.widthfact);
nlen = length(par.length);
periods = kron(par.periods(:)',ones(1,next*nlen));
width = fn_mult(par.periods(:)',par.widthfact(:));
extents = kron(width(:)',ones(1,nlen));
ellips = fn_mult(width(:)',1./par.length(:));
ellips = ellips(:)';

filter_size = pow2(floor(log2(max(ni,nj))) + 1);

% Gabor parameters
norient = length(angles);
nscales = nper*next*nlen;

% coordinates
[Y X] = meshgrid(-(filter_size)/2:(filter_size)/2-1,-(filter_size)/2:(filter_size)/2-1);
FGabor = zeros(filter_size,filter_size,nscale,norient);

fn_progress('Gabor base',nbr_scales*nbr_orient)
i_total = 0;
for i_orient=1:norient
    theta = angles(i_orient);
    U = X*sind(theta) - Y*cosd(theta);
    V = X*cosd(theta) + Y*sind(theta);
    for i_scale=1:nscales
        % Compute filter
        i_total = i_total+1;
        fn_progress(i_total)
        lambda = periods(i_scale);
        sigma = extents(i_scale);
        gamma = ellips(i_scale);
        ATTENUATION = exp(-(U.^2 + (V*gamma).^2) / sigma^2);
        PHASE = exp(2*pi*1i*U/lambda);
        GABOR = ATTENUATION .* PHASE;
        N_F=sum(sum(GABOR.*conj(GABOR)))/filter_size/filter_size;
        Gabor = fftshift(GABOR)/sqrt(N_F);
        Fk = fft2(Gabor);
        FGabor(:,:,i_scale,i_orient) = Fk;
    end
end

%---
function res = test_estimate(data,FGabor,par)

%-- GABOR BASE PARAMETERS --% 

% % Gabor parameters
% nbr_period = length(periods);
% nbr_extent = length(extents);
% nbr_ellips = length(ellips);
% nbr_scales = max([nbr_period nbr_extent nbr_ellips]);
% if ~all([nbr_period nbr_extent nbr_ellips]==nbr_scales)
%     error('number of periods, extents and ellipcities should all be the same or be 1')
% end
% nbr_orient=length(angles);
% M=(nbr_scales*nbr_orient)*filter_size^2; % overcompleteness
% min_angle = angles(1);
% max_angle = angles(end);
% 
% % coordinates
% [fY fX] = meshgrid(1:filter_size,1:filter_size);
% [Y X] = meshgrid(-(filter_size)/2:(filter_size)/2-1,-(filter_size)/2:(filter_size)/2-1);

angles = par.min_angle + ...
    (0:par.nbr_orient-1)/(par.nbr_orient-1)*(par.max_angle-par.min_angle);
nper = length(par.periods); 
next = length(par.widthfact);
nlen = length(par.length);
nbr_orient = length(angles);
nbr_scales = nper*next*nlen;

filter_size = pow2(floor(log2(max(ni,nj))) + 1);

%-- ESTIMATION PARAMETERS --%

seuilPoids=10 ; %weight threshold to consider the weigthed sum as reliable

img1200b=filtx2(data,15,'hm');
tmax=size(img1200b,1) ;
xmax=size(img1200b,2) ;

wdwatt=8 ; %attenuation en x o t: taille de la zone de transition (i.e. du raised cosine)
% wdw=filter_size;    %attenuation 8 pixels bords donnees
att=(.5-.5*cos([1:wdwatt]'/wdwatt*pi));
img1200b(1:wdwatt,:)=img1200b(1:wdwatt,:).*(att*ones(1,xmax)) ;
img1200b(tmax+1-wdwatt:tmax,:)=img1200b(tmax+1-wdwatt:tmax,:).*(att(wdwatt:-1:1)*ones(1,xmax)) ;
szdata=size(img1200b) ;
att=(.5-.5*cos([1:wdwatt]/wdwatt*pi)).^1 ;
img1200b(:,1:wdwatt)=img1200b(:,1:wdwatt).*(ones(tmax,1)*att) ;
img1200b(:,szdata(2)+1-wdwatt:szdata(2))=img1200b(:,szdata(2)+1-wdwatt:szdata(2)).*(ones(tmax,1)*att(wdwatt:-1:1)) ;
myim = zeros(filter_size);
tstart = floor((filter_size - tmax)/2);
tend = tstart + tmax - 1;
xstart = floor((filter_size - xmax)/2);
xend = xstart + xmax - 1;
myim(tstart:tend, xstart:xend) = img1200b;
FI=fft2(myim);
Ic=zeros(tmax,xmax,nbr_scales,nbr_orient);

%-- ESTIMATION LOOP --%

% determine in which quarter we are
If = abs(FI);
take = filter_size/2;
balance = mean(mean(If(2:take,2:take))) - mean(mean(If(2:take,take+2:2*take)));
if balance>0
    % "other" quarter
    disp('detected negative speed')
    angles = angles + 90;
else
    disp('detected positive speed')
end

fn_progress('Gabor weight',nbr_scales*nbr_orient)
i_total = 0;
for i_orient=1:nbr_orient
    theta = angles(i_orient);
    U = X*sind(theta) - Y*cosd(theta);
    V = X*cosd(theta) + Y*sind(theta);
    for i_scale=1:nbr_scales
        % Compute filter
        i_total = i_total+1;
        fn_progress(i_total)
        lambda = periods(i_scale);
        sigma = extents(i_scale);
        gamma = ellips(i_scale);
        ATTENUATION = exp(-(U.^2 + (V*gamma).^2) / sigma^2);
        PHASE = exp(2*pi*1i*U/lambda);
        GABOR = ATTENUATION .* PHASE;
        N_F=sum(sum(GABOR.*conj(GABOR)))/filter_size/filter_size;
        Gabor = fftshift(GABOR)/sqrt(N_F);
        Fk = fft2(Gabor);

        % Estimate
        N_win=512;
        tmpim=ifft2(FI.*conj(Fk))/N_win;
        Ic(:,:,i_scale,i_orient)=tmpim(tstart:tend, xstart:xend);
    end;
end;

coef=abs(Ic).^2 ;
tam=size(coef) ;
tmax=size(coef,1) ;
xmax=size(coef,2) ;

coef0=coef; % save the energies before ponderation
meannrj=mean(mean(mean(coef(1+wdwatt:tmax-wdwatt,1+wdwatt:xmax-wdwatt,:,:),1),2),4) ;
for sc=1:tam(3), %ponderar por la energia media en el canal
    coef(:,:,sc,:)=coef(:,:,sc,:)/meannrj(sc) ;
end
% meanpeaknrj=mean(mean(max(coef(1+wdwatt:tmax-wdwatt,1+wdwatt:xmax-wdwatt,:,:),[],4),1),2) ;
%meanpeaknrj=meanpeaknrj(1:tam(3)-1) ;
%mediannrj=mean(coef,4) ; %median for each position/scale across orientations

coef=coef-1 ; %baseline suppression
coef(coef<0)=0; %baseline suppression

orv = angles;

% compute bin width (as suggested per ivo)
bw(1) = (tand((orv(2)+orv(1))/2) - tand(orv(1)) ) * 2;
for or=2:tam(4)-1
    bw(or) = tand((orv(or+1)+orv(or))/2) - tand((orv(or)+orv(or-1))/2);
end;
bw(tam(4)) = ( tand(orv(tam(4))) - tand( (orv(tam(4))+orv(tam(4)-1)) / 2 ) ) * 2;


poids = zeros(tam(1:2));
anglesc = zeros(tam(1:2));
% speedsc = zeros(tam(1:2));
% poidstg = zeros(tam(1:2));
% angletg = zeros(tam(1:2));
% speedtg = zeros(tam(1:2));
for sc=1:tam(3)
    anglesc_syl{sc} = zeros(tam(1:2)); %#ok<*AGROW>
    %     speedsc_syl{sc} = zeros(tam(1:2));
    poids_syl{sc} =  zeros(tam(1:2));
    %     angletg_syl{sc} = zeros(tam(1:2));
    %     speedtg_syl{sc} = zeros(tam(1:2));
    %     poidstg_syl{sc} =  zeros(tam(1:2));
    for or=1:tam(4)
        anglesc_syl{sc} = anglesc_syl{sc} + coef(:,:,sc,or)*orv(or);
        %         speedsc_syl{sc} = speedsc_syl{sc} + coef(:,:,sc,or)*tand(orv(or));
        poids_syl{sc} = poids_syl{sc} + coef(:,:,sc,or);
        %         angletg_syl{sc} = angletg_syl{sc} + coef(:,:,sc,or)*orv(or)*bw(or);
        %         speedtg_syl{sc} = speedtg_syl{sc} + coef(:,:,sc,or)*tand(orv(or))*bw(or);
        %         poidstg_syl{sc} = poidstg_syl{sc} + coef(:,:,sc,or)*bw(or);
    end
    anglesc = anglesc + anglesc_syl{sc};
    %     speedsc = speedsc + speedsc_syl{sc};
    anglesc_syl{sc} = anglesc_syl{sc}./(poids_syl{sc}+1e-33);
    %     speedsc_syl{sc} = speedsc_syl{sc}./(poids_syl{sc}+1e-33);
    poids = poids + poids_syl{sc};
    
    %     angletg = angletg + angletg_syl{sc};
    %     speedtg = speedtg + speedtg_syl{sc};
    %     angletg_syl{sc} = angletg_syl{sc}./(poidstg_syl{sc}+1e-33);
    %     speedtg_syl{sc} = speedtg_syl{sc}./(poidstg_syl{sc}+1e-33);
    %     poidstg = poidstg + poidstg_syl{sc};
    
end
disp(['Orientation not defined in ',int2str(sum(sum(poids<=0))), ' points on ',int2str(prod(tam(1:2)))]) ; drawnow
disp(['Orientation unreliably defined in ',int2str(sum(sum(poids<seuilPoids))), ' points on ',int2str(prod(tam(1:2)))]) ; drawnow
anglesc=anglesc./(poids+1e-33) ; %weighted sum 
% angletg=angletg./(poidstg+1e-33) ; %weighted sum 
% speedsc=speedsc./(poids+1e-33) ; %weighted sum
% speedtg=speedtg./(poidstg+1e-33) ; %weighted sum 

poids_orig = poids;
anglesc_orig = anglesc;


%%%%%%%%%%%%%%  averaging for the less defined coef.
disp('Propagating neighboring values to pixels for which the angle is undefined or unreliable:'), drawnow
% cvK= [1    10    45   120   210   252   210   120    45    10     1]/1024 ; 
cvK= [1     8    28    56    70    56    28     8     1]/256 ; %[1 6 15 20 15 6 1]/64 ; %[1 4 6 4 1]
it=0;
while it<30 || ~all(poids(:))
    it = it+1;
    anglescM=conv2(anglesc.*poids,cvK'*cvK,'same') ; %convolucion con los vecinos
    poidsM=conv2(poids,cvK'*cvK,'same') ; %convolucion con los vecinos
    anglescM=anglescM./(poidsM+1e-33) ; %recalculo por pesos
    anglesc(poids<seuilPoids)=anglescM(poids<seuilPoids) ;%reemplzar los mal definidos
    poids(poids<seuilPoids)=poidsM(poids<seuilPoids) ; %reemplzar los mal definidos
    maldef=sum(sum(poids<seuilPoids)) ; %number of points unreliably defined
    if maldef==0, break ; %every point is already  accuretely definied
    end
end
disp([int2str(it),' : orientation unreliably defined in ',int2str(maldef), ' points on ',int2str(prod(tam(1:2)))]) ;
    

res.poids = poids;
res.poids_orig = poids_orig;
res.anglesc = anglesc;
res.anglesc_orig = anglesc_orig;
res.poids_syl = poids_syl;
res.anglesc_syl = anglesc_syl;
% res.speedsc = speedsc;
% res.speedsc_syl = speedsc_syl;

% res.poidstg = poidstg;
% res.angletg = angletg;
% res.poidstg_syl = poidstg_syl;
% res.angletg_syl = angletg_syl;
% res.speedtg = speedtg;
% res.speedtg_syl = speedtg_syl;

res.bw = bw; 

% % save Ic and coef
% disp('saving Ic and coef')
% res.Ic = Ic;
% res.coef = coef;

% find for each pixel the best filter
disp('saving best angle'), drawnow
[nt nx nsc ndir] = size(coef0);
% n = nt*nx;
c = coef0;
c = reshape(c,[nt*nx nsc*ndir]);
[dum idx] = max(c,[],2);
idx = reshape(idx,nt,nx);
% isc = 1+mod((idx-1),nsc);
idir = 1+floor((idx-1)/nsc);
angmin = min_angle;
angmax = max_angle;
ang = angmin + (0:(nbr_orient-1))/(nbr_orient-1) * (angmax-angmin);
% ind = 1:nt*nx;
% f = find(idir(:)~=1);
% recons = zeros(nt,nx);
% recons(ind(f)) = Ic(ind(f) + n*(idx(ind(f))-1));
% res.idx = idx;
% res.isc = isc;
% res.idir = idir;
res.bestang = ang(idir);
