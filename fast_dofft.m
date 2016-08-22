function [edges frameidx] = fast_dotensor(edges,params,flag)
% function [edges, frameidx] = fast_dotensor(edges,params,flag)
%---
% Blood flow estimation from fast blood volume measure
% 
% Input:
% - edges   structure with field 'data' containing time courses in edges
% - params  structure ; fields are 
%               .interpt, interpx       resampling
%               .hight, highx           temporal and spatial high pass filtering 
%                                       (to remove heart beat) - done using 
%                                       filty2 and filtx2
%               .lowt, lowx             temporal and spatial smoothing (to enhance lines)
%                                       - done using imfilter
%               .windowwidth            temporal width of the shifting window used for fft2
%               .windowcover            how many shifted window will cover the same poins
%               .sigmaangle             smoothing of the angles empirical distribution
%               .searchangle            half diameter of bla bla (in pi*rad)
%           use fast_params('fft',...) to get the default structure[!not implemented yet!]
% - flag    string ; can be : 
%
% Output:
% - edges   
% - frameidx        indices of frames for which flow has been estimated

warning off MATLAB:divideByZero

% parameters
interpstep = params.interpstep;
hight = params.hight;
highx = params.highx;
lowt = params.lowt;
lowx = params.lowx;
nf = params.windowwidth;
ncover = params.windowcover;
sigma = params.sigmaangle;
searchangle = params.searchangle;

nedges = prod(size(edges));

% go!
fn_progress('flow estimate edge',nedges)
for i=1:nedges
    fn_progress(i)

    y = edges(i).data;
    [nt np] = size(y);

    % high-pass (remove heart beat)
    if hight, y = filty2(y(:,:),hight,'hm'); end
    if highx, y = filtx2(y(:,:),highx,'hm'); end

%     % interpolation
%     [ii jj] = meshgrid(1:interpstep:nx,1:interpstep:nt);
%     y=interp2(y,ii,jj);
% 
%     % low-pass (smoothing to enhance lines)
%     lowt = max(lowt,.1);
%     lowx = max(lowx,.1);
%     h = fspecial('gaussian',[10*ceil(lowt)+1 1],lowt)*fspecial('gaussian',[1 10*ceil(lowx)+1],lowx);
%     y = imfilter(y,h,'replicate');

    % size and number of shifting windows
    dframe = floor(nf/ncover);
    frameidx = floor(nf/2):dframe:floor(nt-nf/2);
    nwindows = length(frameidx);

    % angles in a window
    center = ceil(([nf np]+1)/2);
    [pp ff]=meshgrid((1:np)-center(2),(1:nf)-center(1));
    theta = atan(pp./ff);
    [ang ord]=sort(theta(:)); 
    ord(isnan(ang))=[]; 
    ang(isnan(ang))=[]; 
    reso=1000; 
    angles = pi*(-1/2:1/reso:1/2);
    nangles = length(angles);
    anglesref = 1e-6*ones(1,nangles);
    for j=1:length(ord)
        idx=round(((ang(j)+(pi/2))/pi)*reso)+1; 
        anglesref(idx)=anglesref(idx)+1; 
    end
    
    % loop on windows
    flow = zeros(nwindows,1);
    for k=1:nwindows
        
        % window data
        yw = y((k-1)*dframe+(1:nf),:);
        
        % Fourrier transform
        yf = abs(fft2(yw));
        syf = fftshift(yf);

        % angles distribution
        anglesdist=zeros(1,nangles);
        for j=1:length(ord)
            idx=round(((ang(j)+(pi/2))/pi)*reso)+1; 
            anglesdist(idx)=anglesdist(idx)+syf(ord(j));
        end
        
        % filtering angle distribution
        sanglesdist=filtx(anglesdist,sigma,'lm');
        sanglesref=filtx(anglesref,sigma,'lm');
        tuning = sanglesdist./sanglesref;
    
        % find max and choose the mean value close to this max that best
        % fits the local distribution (tricky)
        [dum argmax] = max(tuning);
        % refine: first choose an interval
        dum = round(reso*searchangle);
        mw = max(1,argmax-dum); Mw = min(nangles,argmax+dum);
        interval = mw:Mw;
        %% reject the parts of the interval where the tuning goes back up
        %tslope = diff(tuning(interval));
        %cw = argmax-mw+1;
        %m=max(find(tslope(1:cw-1)<0)); if isempty(m), m=1; end
        %M=min(find(tslope(cw:end)>0))+cw-1; if isempty(M), M=Mw-mw+1; end
        %interval = interval(m:M);
        %% restrict the interval so that tunings are equal at the extremities
        %if tuning(interval(1))>tuning(interval(end))
        %    interval = interval(1:max(find(tuning(interval)>tuning(interval(1)))));
        %else
        %    interval = interval(min(find(tuning(interval)>tuning(interval(end)))):end);
        %end
        %% then find the mean in this interval
        moyenne = round(sum(interval.*tuning(interval))./sum(tuning(interval)));
        meanangle = angles(moyenne);
        meanangle = angles(argmax);
        
        % deduce the flow
        speed = 1/tan(meanangle)*(np/nf);
        flow(k) = speed;
        
%         figure(3), imagesc(yw), axis image, colormap gray
%         figure(4), imagesc(syf), axis image
%         figure(6), plot(angles,tuning), 
%         pause(.2)
        
    end 
        
    % save
    edges(i).flow = flow;
%         [iii jjj] = meshgrid(1:nx,1:nt);
%         edges(i).dataf = interp2(ii,jj,y,iii,jjj);
%         edges(i).ytyt = interp2(ii,jj,ytyt,iii,jjj);
%         edges(i).ytyx = interp2(ii,jj,ytyx,iii,jjj);
%         edges(i).yxyx = interp2(ii,jj,yxyx,iii,jjj);
    
end

% % optional averaging over same conditions
% if any(strfind(flag,'m'))
%     [nrep ncond] = size(edges);
%     if nrep==1, edges = edges'; [nrep ncond] = size(edges); end
%     fn_progress('averaging conditions',ncond)
%     for k=1:ncond
%         fn_progress(k)
%         E(k).ytyt = fn_means(edges(:,k).ytyt);
%         E(k).ytyx = fn_means(edges(:,k).ytyx);
%         E(k).yxyx = fn_means(edges(:,k).yxyx);
%     end
%     edges = E;
%     nedges = ncond;
% end

warning on MATLAB:divideByZero
