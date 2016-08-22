function edges = fast_dotensor(edges,params,flag)
% function edges = fast_dotensor(edges,params,flag)
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
%               .sigmat, sigmax         temporal and spatial smoothing of
%                                       the structure tensor
%           use fast_params to get the default structure
% - flag    string ; can be : 
%               'A'     only fields dataf, ytyt, ytyx, yxyx are computed
%               'f'     only field flux is computed
%               'A2f'   computes field flux; ytyt, ytyx, yxyx should be
%                       already computed
%               'Af'    computes everything
%               'u'     saves also fu, fv
%               'Am','fm',...   instead of doing the same computations for
%                               every element of the structure, computes a
%                               mean (over columns if it is a 2D structure)
%
% Output:
% - edges   


% parameters
interpstep = params.interpstep;
hight = params.hight; 
highx = params.highx;
% lowt = params.lowt; 
% lowx = params.lowx; 
% sigmat = params.sigmat; 
% sigmax = params.sigmax; 
% hight = params.hight/interpstep;
% highx = params.highx/interpstep;
lowt = params.lowt/interpstep;
lowx = params.lowx/interpstep;
sigmat = params.sigmat/interpstep;
sigmax = params.sigmax/interpstep;

nedges = prod(size(edges));

% compute the structure tensor
if ~any(strfind(flag,'2'))
    fn_progress('preprocessing',nedges)
    for i=1:nedges
        fn_progress(i)
        
        y = edges(i).data;  
        [nt nx] = size(y);
        
        % high-pass (remove heart beat)
        if hight, y = filty2(y(:,:),hight,'hm'); end
        if highx, y = filtx2(y(:,:),highx,'hm'); end
        
% ESSAIS INFRUCTUEUX
%         % interpolation
%         [ii jj] = meshgrid(1:interpstep:nx,1:interpstep:nt);
%         y=interp2(y,ii,jj);
%         
%         % interpolation of a dirac
%         [ii jj] = meshgrid(1:interpstep:ceil(10*lowt),1:interpstep:ceil(10*lowt));
%         [iii jjj] = meshgrid(1:interpstep:ceil(10*lowt),1:interpstep:ceil(10*lowt));
%         d = zeros(ceil(10*lowt)); d(ceil(5*lowt),ceil(5*lowt))=1;
%         d = interp2(d,ii,jj);
%         
%         % low-pass (smoothing to enhance lines)
%         lowt = max(lowt,.1);
%         lowx = max(lowx,.1);
%         h = fspecial('gaussian',[ceil(6*lowt) 1],lowt)*fspecial('gaussian',[1 ceil(6*lowx)],lowx);
%         hd = conv2(d,h,'same');
%         hd = interp2(ii,jj,hd,iii,jjj);
%         y = imfilter(y,hd,'replicate');
        
        % interpolation
        [ii jj] = meshgrid(1:interpstep:nx,1:interpstep:nt);
        y=interp2(y,ii,jj);

        % low-pass (smoothing to enhance lines)
        lowt = max(lowt,.1);
        lowx = max(lowx,.1);
        h = fspecial('gaussian',[10*ceil(lowt)+1 1],lowt)*fspecial('gaussian',[1 10*ceil(lowx)+1],lowx);
        y = imfilter(y,h,'replicate');
        

        % derivatives
        [yx yt] = gradient(y); 
        ytyt = yt.*yt;
        ytyx = yt.*yx;
        yxyx = yx.*yx;
        
        % average the structure tensor
        h = fspecial('gaussian',[2*ceil(sigmat)+1 1],sigmat)*fspecial('gaussian',[1 2*ceil(sigmax)+1],sigmax);
        ytyt = imfilter(ytyt,h,'replicate');
        ytyx = imfilter(ytyx,h,'replicate');
        yxyx = imfilter(yxyx,h,'replicate');
        
        % save
        if any(strfind(flag,'A')) && ~any(strfind(flag,'A2f'))
            [iii jjj] = meshgrid(1:nx,1:nt);
            edges(i).dataf = interp2(ii,jj,y,iii,jjj);
            edges(i).ytyt = interp2(ii,jj,ytyt,iii,jjj);
            edges(i).ytyx = interp2(ii,jj,ytyx,iii,jjj);
            edges(i).yxyx = interp2(ii,jj,yxyx,iii,jjj);
        end
        
    end  
end

% optional averaging over same conditions
if any(strfind(flag,'m'))
    [nrep ncond] = size(edges);
    if nrep==1, edges = edges'; [nrep ncond] = size(edges); end
    fn_progress('averaging conditions',ncond)
    for k=1:ncond
        fn_progress(k)
        E(k).ytyt = fn_means(edges(:,k).ytyt);
        E(k).ytyx = fn_means(edges(:,k).ytyx);
        E(k).yxyx = fn_means(edges(:,k).yxyx);
    end
    edges = E;
    nedges = ncond;
end

% estimate the flow
if any(strfind(flag,'f'))
    fn_progress('compute flow',nedges)
    warning off MATLAB:divideByZero
    for i=1:nedges
        fn_progress(i)
        
        ytyt = edges(i).ytyt;
        ytyx = edges(i).ytyx;
        yxyx = edges(i).yxyx;
        [nt nx] = size(ytyt);
        
        % compute its eigenvectors
        fu = zeros(nt,nx); fv = zeros(nt,nx); %poids = zeros(nt,np);
        for k=1:nt*nx
            A = [ytyt(k) ytyx(k); ytyx(k) yxyx(k)];
            [u s v] = svd(A);
            poids = s(1,1)/s(2,2);
            fu(k) = u(1,1); fv(k) = u(2,1);
            fv(k) = fv(k);
        end
        
        % % if one wants to save fu and fv
        % norme = sqrt(fu.^2 + fv.^2);
        % edges(i).fu = fu./norme;
        % edges(i).fv = fv./norme;
        
        % compute the flow
        if any(strfind(flag,'u'))
            edges(i).fu = fu;
            edges(i).fv = fv;
        end
        edges(i).flux = -(fu./fv);
    end    
    warning on MATLAB:divideByZero
end
