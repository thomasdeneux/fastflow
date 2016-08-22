function [Y pos] = fast_syl_NEWFAKE_tom(thresh,step,jitter,shotnoise,ni,nj,nt)
% function [Y pos] = fast_syl_NEWFAKE_tom(thresh,step,jitter,shotnoise,ni,nj,nt)
%---
%
% this function creates a simulation of data Y. Red Blood Cells patterns are simulated
% as gaussians that move with speed STEP pixels/frame until frame begchange-1, 
% then switch to step+change pix/frame until frame 600, and then go back to step pix/frame. 
% The positions of the gaussians are jittered by JITTER, which is in turn scaled by the speed (step).
% THRESH (betwen 0 and 1) controls the density of RBCs (the higher, the lower the
% density). 
% All input values can be fractional
% 
% Output:
% - Y       nj*ni*nt array: movie, the vessel is a column of a frame, i.e.
%           it has nj points
% - pos     nrbcs*nt array: positions of the RBCs
% 
% Ivo, 23 June 2006
% Thomas, 21 July 2009, 14 September 2010


%step=1;
change=step/20; 
%change = 0.1;
%thresh=0.9;
gaussamp=20; %amplitude (contrast) of RBCs : the real stuff
%gaussamp=200; %amplitude (contrast) of RBCs: to make life easy
gausslat=5; % x-dimension of RBCs (total x-size will be 2*gausslat+1, the y-size will be determined by an analytical gaussian expression)
sigma=1;  % lateral fall-off of RBCs (i.e. sigma: RBCs are modelled as a gaussian)

begchange=101;
endchange=200;
Y0=zeros(nj,ni,nt);
% creation of "shot-noise" in pseudo-data images Y
for i=1:nt
    randn('state',i+sum(100*clock));
    Y0(:,:,i)=shotnoise*randn(nj,ni);
end
%Y0=Y; %the non-jittered one
%Y1=Y; %the uniformly-jittered one
% creating the shape of a RBC
h=gaussamp*fspecial('gaussian',[1,2*gausslat+1],sigma);

% generating a random distribution of RBC positions (a vector with values
% between -nreserve and +nj)
nreserve = ceil((nt-1)*(step+change));
a = rand(1,nj+nreserve);
indices = find(a>thresh) - nreserve;
nRBC = length(indices);

% generating the temporal evolution

Y = Y0;
pos = zeros(length(indices),nt+1);
pos(:,1) = indices+step;

for t = 1:nt %loop over time
    % update position of the RBCs
    % (regular motion)
    if t<begchange || t>endchange
        indices = indices + step;
    else
        indices = indices + (step+change);
    end
    % (jitter)
    indices = indices + jitter*step*randn(1,nRBC);
     
    % memorize the positions
    pos(:,t+1) = indices;

    % 'print' each RBC in the movie frame
    vec = zeros(nj,2*gausslat+1);       
    for i = 1:nRBC
        vec = vec + exp(-(((1:nj)'-indices(i))/sigma).^2)*h;
    end
    Y(:,round(ni/2)-gausslat:round(ni/2)+gausslat,t)=  Y(:,round(ni/2)-gausslat:round(ni/2)+gausslat,t)-vec;
end


Y=Y+100;


