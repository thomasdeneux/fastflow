function fake = fast_FAKE_thlim(thresh,step,jitter,ni,nj,nt);

% [fakega,fakeno] = NEWFAKEanlwithrand4TEST_syl(thresh,step,jitter,ni,nj,nt);
%
% this function creates a simulation of data Y. Red Blood Cells patterns are simulated
% as gaussians that move with speed STEP pixels/frame until frame begchange-1, 
% then switch to step+change pix/frame until frame 600, and then go back to step pix/frame. 
% The positions of the gaussians are jittered by JITTER, which is in turn scaled by the speed (step).
% THRESH (betwen 0 and 1) controls the density of RBCs (the higher, the lower the
% density). 
% All input values can be fractional
% 
% Ivo, 23 June 2006

%step=1;
change=step/10; 
%change = 0.1;
%thresh=0.9;
%gaussamp=20; %amplitude (contrast) of RBCs : the real stuff
gaussamp=200; %amplitude (contrast) of RBCs: to make life easy
gausslat=5; % x-dimension of RBCs (total x-size will be 2*gausslat+1, the y-size will be determined by an analytical gaussian expression)
sigma=1;  % lateral fall-off of RBCs (i.e. sigma: RBCs are modelled as a gaussian)

% those are in parameters now!!! (Sylvain Takerkart)
nj=100; % y-dimension: the RBCs shall move in a vessel oriented vertically
ni=100; % x-dimension
nt=300; % number of frames



%jitter=1/3;
begchange=101;
endchange=200;
Y0=zeros(nj,ni,nt);
% creation of "shot-noise" in pseudo-data images Y
for i=1:nt
    randn('state',i+sum(100*clock));
    Y0(:,:,i)=randn(nj,ni);
end
%Y0=Y; %the non-jittered one
%Y1=Y; %the uniformly-jittered one
% creating the shape of a RBC
h=gaussamp*fspecial('gaussian',[1,2*gausslat+1],sigma);
% generating a random distribution of RBC positions (a vertical vector)
indices=[];
while (size(indices,1) == 0)
    rand('state',sum(100*clock));
    a = rand(nj+ceil((nt-1)*(step+change)),1);
    a=a>thresh;
    [indices,dummy]=find(a);
end

% generating the temporal evolution

for j=1:length(jitter)
    Y = Y0;
    indices1=indices;
    speedband=zeros(nj,2*gausslat+1,nt-1);
    for t = 1:nt %loop over time
        
        rr=jitter(j)*step*randn(size(indices1));
        indices1=indices1+rr;
        vec = zeros(nj,2*gausslat+1);
        speed=zeros(nj,2*gausslat+1);
        for i = 1:max(size(indices)) % loop over all RBCs
            if (t<begchange)
                vec(:,:)=vec(:,:) + exp( -( ( [1:nj]' + (t-1)*(step) - indices1(i) )/sigma  ).^2  )*h;
                vitesse = step+rr(i);
                speed(:,:)=speed(:,:) + exp( -( ( [1:nj]' + (t-1)*(step) - indices1(i) )/sigma  ).^2  )*h*vitesse;
            elseif (t>begchange-1)&(t<endchange+1)
                vec(:,:)=vec(:,:)+exp(-(([1:nj]'+((begchange-2)*step+(t-(begchange-1))*(step+change))-indices1(i))/sigma).^2)*h;
                vitesse = step+rr(i);
                speed(:,:)=speed(:,:) + exp(-(([1:nj]'+((begchange-2)*step+(t-(begchange-1))*(step+change))-indices1(i))/sigma).^2)*h*vitesse;
            else
                vec(:,:)=vec(:,:)+exp(-(([1:nj]'+((begchange-2)*step+(endchange-begchange)*(step+change)+(t-(endchange-1))*step)-indices1(i))/sigma).^2)*h;
                vitesse = step+rr(i);
                speed(:,:)=speed(:,:) + exp(-(([1:nj]'+((begchange-2)*step+(endchange-begchange)*(step+change)+(t-(endchange-1))*step)-indices1(i))/sigma).^2)*h*vitesse;           
            end

        end
        speedband(:,:,t)=speed/sum(sum(vec)); 
        Y(:,round(ni/2)-gausslat:round(ni/2)+gausslat,t)=  Y(:,round(ni/2)-gausslat:round(ni/2)+gausslat,t)-vec;
    end
    Y=Y+100;
    %Y0=Y0+100;
    
    CS=zeros(nj,ni);
    CS(:,round(ni/2)-gausslat:round(ni/2)+gausslat)=CS(:,round(ni/2)-gausslat:round(ni/2)+gausslat)-h(ones(nj,1),:);
    [mask CSU CSV] = fast_structure(CS+1000+randn(100,100),'csu');
    
    fake{j}.Y = Y;
    fake{j}.CS = CS;
    fake{j}.CSU = CSU;
    fake{j}.CSV=CSV;
    fake{j}.speed=speedband;

end;


