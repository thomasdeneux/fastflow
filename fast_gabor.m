function [x s] = fast_gabor(varargin)
% [x s] = fast_gabor(data[,par[,filterbase]])
% filterbase = fast_gabor([par,]ni,nj)
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
% Thomas Deneux, 2008, 2009, 2010

if nargin==1 && isequal(varargin{1},'par')
    % default parameters
    x = defaultparameters;
elseif nargin==2
    % build base of filters
    par = defaultparameters;
    [ni nj] = deal(varargin{:});
    x = test_base(par,ni,nj);
elseif nargin==3 && isstruct(varargin{1})
    % build base of filters
    [par ni nj] = deal(varargin{:});
    x = test_base(par,ni,nj);
else
    % estimate 
    data = varargin{1};
    if nargin>=2
        par = varargin{2}; 
    else
        par = defaultparameters; 
    end
    if nargin>=3
        FGabor = varargin{3};
    else
        [ni nj] = size(data);
        FGabor = test_base(par,ni,nj);
    end
    [x s] = test_estimate(data,FGabor,par);
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
% par.gabordisplay = false;

%---
function FGabor = test_base(par,ni,nj)

persistent basesave parsave nisave njsave 
if isempty(nisave), nisave=0; njsave=0; end

% do not recompute filter if parameters haven't changed
tmp = fn_structmerge(defaultparameters,par,'skip','type');
if isequal(tmp,parsave) && ni<=nisave && nj<=njsave
    FGabor = basesave;
    return
end
parsave = par;
nisave  = max(ni,nisave);
njsave  = max(nj,njsave);

% some more parameters
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

% update ni, nj
minsize = 2*ceil(2*max([width(:)' par.length(:)']))+1;
ni = max(nisave,minsize);
nj = max(njsave,minsize);

% Gabor parameters
norient = length(angles);
nscales = nper*next*nlen;

% coordinates
[Y X] = meshgrid(-nj/2:nj/2-1,-ni/2:ni/2-1);
FGabor = zeros(ni,nj,nscales,norient);
% display
% if par.gabordisplay
nshow = 5;
show_orient = false(1,norient);
show_orient(round(linspace(1,norient,nshow))) = true;
show_size = ceil(2*max([width(:)' par.length(:)']));
show_i = ceil(ni/2) + (-show_size:show_size);
show_j = ceil(nj/2) + (-show_size:show_size);
show_n = 2*show_size+1;
ShowGabor = zeros(show_n,show_n,nscales,nshow);
% end

fn_progress('Gabor base',nscales*norient)
i_total = 0;
for i_orient=1:norient
    theta = angles(i_orient);
    U = X*sind(theta) - Y*cosd(theta);
    V = X*cosd(theta) + Y*sind(theta);
    for i_scale=1:nscales
        % Compute filter
        i_total = i_total+1;
        if ~mod(i_total,10), fn_progress(i_total), end
        lambda = periods(i_scale);
        sigma = extents(i_scale);
        gamma = ellips(i_scale);
        ATTENUATION = exp(-(U.^2 + (V*gamma).^2) / sigma^2);
        PHASE = exp(2*pi*1i*U/lambda);
        GABOR = ATTENUATION .* PHASE;
        N_F=sum(sum(GABOR.*conj(GABOR)))/ni/nj;
        Gabor = fftshift(GABOR)/sqrt(N_F);
        if show_orient(i_orient)
            ShowGabor(:,:,i_scale,sum(show_orient(1:i_orient))) = ...
                real(GABOR(show_i,show_j));
        end
        Fk = fft2(Gabor);
        FGabor(:,:,i_scale,i_orient) = Fk;
    end
end

% informative display
% if par.gabordisplay
ShowGabor = reshape(permute(ShowGabor,[1 4 2 3]),show_n*nshow,show_n*nscales);
if ishandle(940)
    ha = findobj(940,'type','axes');
else
    figure(940)
    set(940,'numbertitle','off','name','fast_gabor: filter set (all scales, selected orientations)')
    ha = axes('position',[0 0 1 1]);
    colormap gray
end
imagesc(ShowGabor,'parent',ha)
% end

% save base in persistent memory
basesave = FGabor;

%---
function [x res] = test_estimate(data,FGabor,par)

angles = par.min_angle + ...
    (0:par.nbr_orient-1)/(par.nbr_orient-1)*(par.max_angle-par.min_angle);
nper = length(par.periods); 
next = length(par.widthfact);
nlen = length(par.length);
nbr_orient = length(angles);
nbr_scales = nper*next*nlen;

[ni0 nj0] = size(data);
[ni nj dum1 dum2] = size(FGabor); %#ok<NASGU>
istart = 1 + floor((ni-ni0)/2);
iend = ni0 + floor((ni-ni0)/2);
jstart = 1 + floor((nj-nj0)/2);
jend = nj0 + floor((nj-nj0)/2);

seuilPoids=10 ; %weight threshold to consider the weigthed sum as reliable

myim = ones(ni,nj)*mean(data(:)); % outside value equal to mean -> minimize edge effects
myim(istart:iend,jstart:jend) = data;
FI=fft2(myim);
Ic=zeros(ni0,nj0,nbr_scales,nbr_orient);

% determine in which quarter we are
If = abs(FI);
take = floor([ni nj]/2);
balance = mean(mean(If(2:take(1),2:take(2)))) - mean(mean(If(2:take(1),nj-take(2)+2:nj)));
if balance>0
    % "other" quarter
    disp('detected negative speed')
    angles = - angles;
    myim = flipud(myim);
    FI=fft2(myim);
else
    disp('detected positive speed')
end

%-- ESTIMATION LOOP --%

fn_progress('Gabor weight',nbr_scales*nbr_orient)
i_total = 0;
for i_orient=1:nbr_orient
    for i_scale=1:nbr_scales
        % filter
        i_total = i_total+1;
        if ~mod(i_total,10), fn_progress(i_total), end
        Fk = FGabor(:,:,i_scale,i_orient);

        % Estimate
        N_win=512;
        tmpim=ifft2(FI.*conj(Fk))/N_win;
        Ic(:,:,i_scale,i_orient)=tmpim(istart:iend,jstart:jend);
    end;
end;

coef=abs(Ic).^2 ;
tam=size(coef) ;

coef0=coef; % save the energies before ponderation
meannrj=mean(mean(mean(coef,1),2),4) ;
for sc=1:tam(3), %ponderar por la energia media en el canal
    coef(:,:,sc,:)=coef(:,:,sc,:)/meannrj(sc) ;
end

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
for sc=1:tam(3)
    anglesc_syl{sc} = zeros(tam(1:2)); %#ok<*AGROW>
    poids_syl{sc} =  zeros(tam(1:2));
    for or=1:tam(4)
        anglesc_syl{sc} = anglesc_syl{sc} + coef(:,:,sc,or)*orv(or);
        poids_syl{sc} = poids_syl{sc} + coef(:,:,sc,or);
    end
    anglesc = anglesc + anglesc_syl{sc};
    anglesc_syl{sc} = anglesc_syl{sc}./(poids_syl{sc}+1e-33);
    poids = poids + poids_syl{sc};
end
% disp(['Orientation not defined in ',int2str(sum(sum(poids<=0))), ' points on ',int2str(prod(tam(1:2)))]) ; drawnow
% disp(['Orientation unreliably defined in ',int2str(sum(sum(poids<seuilPoids))), ' points on ',int2str(prod(tam(1:2)))]) ; drawnow
maldef = sum(poids(:)<seuilPoids);
npoints = prod(tam(1:2));
fprintf('reliability: %.0f%% initially -> \n',100-100*maldef/npoints), drawnow
anglesc=anglesc./(poids+1e-33) ; %weighted sum 

poids_orig = poids;
anglesc_orig = anglesc;


%%%%%%%%%%%%%%  averaging for the less defined coef.
% disp('Propagating neighboring values to pixels for which the angle is undefined or unreliable:'), drawnow
cvK= [1     8    28    56    70    56    28     8     1]/256 ; 
it=0;
while it<30 || ~all(poids(:))
    it = it+1;
    anglescM=conv2(anglesc.*poids,cvK'*cvK,'same') ; %convolucion con los vecinos
    poidsM=conv2(poids,cvK'*cvK,'same') ; %convolucion con los vecinos
    anglescM=anglescM./(poidsM+1e-33) ; %recalculo por pesos
    anglesc(poids<seuilPoids)=anglescM(poids<seuilPoids) ;%reemplzar los mal definidos
    poids(poids<seuilPoids)=poidsM(poids<seuilPoids) ; %reemplzar los mal definidos
    maldef=sum(poids(:)<seuilPoids) ; %number of points unreliably defined
    if maldef==0, break ; %every point is already  accuretely definied
    end
end
%disp([int2str(it),' : orientation unreliably defined in ',int2str(maldef), ' points on ',int2str(prod(tam(1:2)))]) ;
fprintf('\b%.0f%% after %i propagation iterations\n',100-100*maldef/npoints,it), drawnow

% re-reverse the results if negative speeds
if balance>0
    poids = flipud(poids);
    poids_orig = flipud(poids_orig);
    anglesc = flipud(anglesc);
    for i=1:nbr_scales, poids_syl{i} = flipud(poids_syl{i}); end
    for i=1:nbr_scales, anglesc_syl{i} = flipud(anglesc_syl{i}); end
    bw = flipud(bw);
end

res.poids = poids;
res.poids_orig = poids_orig;
res.anglesc = anglesc;
res.anglesc_orig = anglesc_orig;
res.poids_syl = poids_syl;
res.anglesc_syl = anglesc_syl;
res.bw = bw; 

% find for each pixel the best filter
[nt nx nsc ndir] = size(coef0);
c = coef0;
c = reshape(c,[nt*nx nsc*ndir]);
[dum idx] = max(c,[],2);
idx = reshape(idx,nt,nx);
idir = 1+floor((idx-1)/nsc);
bestang = angles(idir);
if balance>0
    bestang = flipud(bestang);
end
res.bestang = bestang;

% velocity result
x = cotd(res.anglesc);


