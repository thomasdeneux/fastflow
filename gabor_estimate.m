function res = fast_gabor(data,aff,filter_struct)
% [anglesc,angl,coef]=fast_gabor(data,aff,filter_struct)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%         parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seuilPoids=10 ; %weight threshold to consider the weigthed sum as reliable

img1200b=filtx2(data,15,'hm');
tmax=size(img1200b,1) ;
xmax=size(img1200b,2) ;

if exist('aff')==0, %exist('img1200')==0, %1==1, %charger donnees
    aff=1;
end
% $$$ if exist('elim90')==0, %exist('img1200')==0, %1==1, %charger donnees
% $$$     elim90=1 ;    
% $$$ end
if exist('filter_struct')==0, %1==1, %charger donnees
% $$$     if xmax<64-15, %si les donn�es ont moins de 64 pixels en x
% $$$        fileGabor='GaborV6_px256_sc4_or18.mat'; 
% $$$     elseif xmax<128-15, %si les donn�es ont moins de  128 pixels en x
% $$$        fileGabor='GaborV6_px256_sc4_or18.mat'; 
% $$$     elseif xmax<=256, %si les donn�es ont  moins de 256 pixels en x
% $$$        fileGabor='GaborV6_px256_sc4_or18.mat'; 
% $$$     else
% $$$         disp('Error: x data too long') ;
% $$$     end
% $$$     disp(['loading Gabor filters from: ',fileGabor]) ;
% $$$     eval(['load ',fileGabor, ' ;']) ;
  % load Gaborsc2-5or36
    nbr_scales = 4;
    nbr_orient = 18;
    filter_size = pow2(floor(log2(max(size(img1200b)))) + 1);
    filter_struct = log_gabor_pyramid_syl(filter_size,nbr_scales,nbr_orient);
else
    filter_size = filter_struct.s;
    nbr_scales = filter_struct.n_scale;
    nbr_orient = filter_struct.n_orient;
end 

wdwatt=8 ; %attenuation en x o t: taille de la zone de transition (i.e. du raised cosine)
wdw=filter_struct.s;    %attenuation 8 pixels bords donnees
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
for i_scale=1:nbr_scales
    for i_orient=1:nbr_orient
        %%f_win=ifft2(filter_struct.Fk(:,:,i_scale,i_orient));
        %f_win=filter_struct.Gabor(:,:,i_scale,i_orient);
        %N_win=sqrt(sum(f_win(:).*conj(f_win(:))));
        N_win=512;
        tmpim=ifft2(FI.*conj(filter_struct.Fk(:,:,i_scale,i_orient)))/N_win;
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
meanpeaknrj=mean(mean(max(coef(1+wdwatt:tmax-wdwatt,1+wdwatt:xmax-wdwatt,:,:),[],4),1),2) ;
%meanpeaknrj=meanpeaknrj(1:tam(3)-1) ;
%mediannrj=mean(coef,4) ; %median for each position/scale across orientations

coef=coef-1 ; %baseline suppression
coef(coef<0)=0; %baseline suppression
% determine the band where the energy is max
[tmp middle_band] = max(squeeze(mean(mean(mean(coef,1),2),3)));



for or=1:tam(4)
    orv(or) = (filter_struct.min_angle + (filter_struct.max_angle - filter_struct.min_angle)*(or-1)/(tam(4)-1)) ...
          * 180 / pi;
end;

% compute bin width (as suggested per ivo)
bw(1) = (tand((orv(2)+orv(1))/2) - tand(orv(1)) ) * 2;
for or=2:tam(4)-1
    bw(or) = tand((orv(or+1)+orv(or))/2) - tand((orv(or)+orv(or-1))/2);
end;
bw(tam(4)) = ( tand(orv(tam(4))) - tand( (orv(tam(4))+orv(tam(4)-1)) / 2 ) ) * 2;


poids = zeros(tam(1:2));
anglesc = zeros(tam(1:2));
speedsc = zeros(tam(1:2));
poidstg = zeros(tam(1:2));
angletg = zeros(tam(1:2));
speedtg = zeros(tam(1:2));
for sc=1:tam(3)
    anglesc_syl{sc} = zeros(tam(1:2)); %#ok<*AGROW>
    speedsc_syl{sc} = zeros(tam(1:2));
    poids_syl{sc} =  zeros(tam(1:2));
    angletg_syl{sc} = zeros(tam(1:2));
    speedtg_syl{sc} = zeros(tam(1:2));
    poidstg_syl{sc} =  zeros(tam(1:2));
    for or=1:tam(4)
% $$$         orv = (filter_struct.min_angle + (filter_struct.max_angle - filter_struct.min_angle)*(or-1)/(tam(4)-1)) ...
% $$$               * 180 / pi;
        
        
% $$$         shifted_or = or + floor(tam(4)/2) - middle_band;
% $$$         if shifted_or < 1
% $$$             orv=180/tam(4)*(or-1+.5*mod(sc+1,2)) + 180;
% $$$         elseif shifted_or > tam(4)
% $$$             orv=180/tam(4)*(or-1+.5*mod(sc+1,2)) - 180;
% $$$         else
% $$$             orv=180/tam(4)*(or-1+.5*mod(sc+1,2));
% $$$         end;
        
% $$$         orv=180/tam(4)*(or-1+.5*mod(sc+1,2)) ; %due to the scale shift
% $$$         if (elim90==90)&(or==tam(4)),
% $$$             orv=orv-180 ; %considerar la ultima orientacion como la misma -180�.
% $$$         end
% $$$         if (elim90==180)&(or==1),
% $$$             orv=orv+180 ; %considerar la primera orientacion como la misma+180�.
% $$$         end

        anglesc_syl{sc} = anglesc_syl{sc} + coef(:,:,sc,or)*orv(or);
        speedsc_syl{sc} = speedsc_syl{sc} + coef(:,:,sc,or)*tand(orv(or));
        poids_syl{sc} = poids_syl{sc} + coef(:,:,sc,or);
        angletg_syl{sc} = angletg_syl{sc} + coef(:,:,sc,or)*orv(or)*bw(or);
        speedtg_syl{sc} = speedtg_syl{sc} + coef(:,:,sc,or)*tand(orv(or))*bw(or);
        poidstg_syl{sc} = poidstg_syl{sc} + coef(:,:,sc,or)*bw(or);
    end
    anglesc = anglesc + anglesc_syl{sc};
    speedsc = speedsc + speedsc_syl{sc};
    anglesc_syl{sc} = anglesc_syl{sc}./(poids_syl{sc}+1e-33);
    speedsc_syl{sc} = speedsc_syl{sc}./(poids_syl{sc}+1e-33);
    poids = poids + poids_syl{sc};
    
    angletg = angletg + angletg_syl{sc};
    speedtg = speedtg + speedtg_syl{sc};
    angletg_syl{sc} = angletg_syl{sc}./(poidstg_syl{sc}+1e-33);
    speedtg_syl{sc} = speedtg_syl{sc}./(poidstg_syl{sc}+1e-33);
    poidstg = poidstg + poidstg_syl{sc};
    
end
disp(['Orientation not defined in ',int2str(sum(sum(poids<=0))), ' points on ',int2str(prod(tam(1:2)))]) ;
    disp(['Orientation unreliably defined in ',int2str(sum(sum(poids<seuilPoids))), ' points on ',int2str(prod(tam(1:2)))]) ;  
anglesc=anglesc./(poids+1e-33) ; %weighted sum 
angletg=angletg./(poidstg+1e-33) ; %weighted sum 
speedsc=speedsc./(poids+1e-33) ; %weighted sum 
speedtg=speedtg./(poidstg+1e-33) ; %weighted sum 

poids_orig = poids;
anglesc_orig = anglesc;


%%%%%%%%%%%%%%  averaging for the less defined coef.
disp('Propagating neighboring values to pixels for which the angle is undefined or unreliable:')
cvK= [1    10    45   120   210   252   210   120    45    10     1]/1024 ; 
cvK= [1     8    28    56    70    56    28     8     1]/256 ; %[1 6 15 20 15 6 1]/64 ; %[1 4 6 4 1]
for it=1:30,
    anglescM=conv2(anglesc.*poids,cvK'*cvK,'same') ; %convolucion con los vecinos
    poidsM=conv2(poids,cvK'*cvK,'same') ; %convolucion con los vecinos
    anglescM=anglescM./(poidsM+1e-33) ; %recalculo por pesos
    anglesc(poids<seuilPoids)=anglescM(poids<seuilPoids) ;%reemplzar los mal definidos
    poids(poids<seuilPoids)=poidsM(poids<seuilPoids) ; %reemplzar los mal definidos
    maldef=sum(sum(poids<seuilPoids)) ; %number of points unreliably defined
    if maldef==0, break ; %every point is already  accuretely definied
    end
end
if min(min(poids))==0,
    disp([int2str(it),' : orientation not defined in ',int2str(sum(sum(poids<=0))), ' points on ',int2str(prod(tam(1:2)))]) ;
end
disp([int2str(it),' : orientation unreliably defined in ',int2str(maldef), ' points on ',int2str(prod(tam(1:2)))]) ;
    

res.poids = poids;
res.poids_orig = poids_orig;
res.anglesc = anglesc;
res.anglesc_orig = anglesc_orig;
res.poids_syl = poids_syl;
res.anglesc_syl = anglesc_syl;
res.speedsc = speedsc;
res.speedsc_syl = speedsc_syl;

res.poidstg = poidstg;
res.angletg = angletg;
res.poidstg_syl = poidstg_syl;
res.angletg_syl = angletg_syl;
res.speedtg = speedtg;
res.speedtg_syl = speedtg_syl;

res.bw = bw; 

% % save Ic and coef
% disp('saving Ic and coef')
% res.Ic = Ic;
% res.coef = coef;

% find for each pixel the best filter
disp('saving best angle')
[nt nx nsc ndir] = size(coef0);
n = nt*nx;
c = coef0;
c = reshape(c,[nt*nx nsc*ndir]);
[dum idx] = max(c,[],2);
idx = reshape(idx,nt,nx);
% isc = 1+mod((idx-1),nsc);
idir = 1+floor((idx-1)/nsc);
angmin = filter_struct.min_angle;
angmax = filter_struct.max_angle;
nbr_orient = filter_struct.n_orient;
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



