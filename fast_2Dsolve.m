function [edges, FLUX] = fast_2Dsolve(Y,edges)
% obsolete: use fast_interpalongedges and fast_solve2D instead

% function [edges FLUX] = fast_2Dsolve(Y)
% function [edges FLUX] = fast_2Dsolve(Yf,edges)
%---
% Blood flow estimation from fast blood volume measure
% 
% Input:
% - Y       blood volume (space x space x time)
% - Yf      the same, after spatial masking and smoothing
% - edges   structure with field 'points' containing segmentation of vessels
%           it is estimated if not specified           
%
% Output:
% - edges   two fields are added: 'data' (interpolated data)
%           and 'fu', 'fv' (estimated "directions" in 'data')
% - FLUX    blood flow (time x space x space x 2 sparse array)

error('obsolete: use fast_interpalongedges and fast_solve2D instead')

if nargin<2
    [CSU CSV mask edges lab] = fast_structure(Y(:,:,1));
    Yn = Y ./ repmat(mean(Y,3),[1 1 nt]); clear Y
    Yf = fast_anisosmooth(Yn,CSU,CSV,mask); clear Yn
else 
    Yf = Y; clear Y
    mask = (Yf(:,:,1)~=0);
end



Yf = permute(Yf,[3 1 2]);
    
% data
edges(1).data=[];
disp('interpolating edges          ')
for i=1:nedges
    %fprintf('\b\b\b\b\b\b\b\b\b\b%4i/%4i\n',i,nedges)
    x = edges(i).points2;
    nx = size(x,1);
    data = zeros(nt,nx);
    for k=1:nx
        qi = floor(x(k,1)); ri = x(k,1)-qi;
        qj = floor(x(k,2)); rj = x(k,2)-qj;
        f = find(mask(qj:qj+1,qi:qi+1));
        yf0 = Yf(:,qj:qj+1,qi:qi+1);
        r = [(1-rj) rj]'*[(1-ri) ri]; 
        if sum(r(f))==0
            r(:)=nan;
        else
            r = r/sum(r(f));
        end
        for z=f'
            data(:,k) = data(:,k) + r(z)*yf0(:,z);
        end
    end
    % probl�me des trous : pas r�solu compl�tement
    f = find(isnan(data(1,:)));
    while ~isempty(f)
        for k=f
            kk = [k-1 k+1];
            if k==1, kk=k+1; elseif k==nx, kk=k-1; end
            data(:,k) = nanmean(data(:,kk)')';
        end
        f = find(isnan(data(1,:)));
    end
    edges(i).data = data;
end

% flow smoothing (more exactly: represents the extent of neighbourghood
% when making structure matrice, see below)
h = fspecial('gaussian',[10 10],3);

disp('estimating flow           ')
if nargout>1, FLUX = zeros(nt,nj*ni,1); end
warning off MATLAB:divideByZero
for i=1:nedges
    fprintf('\b\b\b\b\b\b\b\b\b\b%4i/%4i\n',i,nedges)
    yf = edges(i).data;
    yff = imfilter(filty2(yf(:,:),15,'hm'),fspecial('gaussian',[5 1],1)*fspecial('gaussian',[1 5],1));
    
    np = size(yff,2);

%     twidth = 1;
%     xwidth = 6;
%     nxplus = 2*twidth*xwidth+2;
%     yff = [repmat(yff(:,1),1,nxplus) yff repmat(yff(:,end),1,nxplus)]; 
%     flux = zeros(nt,np+2*nxplus);
%     %fprintf('estimation           ')
%     for t=1+twidth:nt-twidth
%         %fprintf('\b\b\b\b\b\b\b\b\b\b%0.4i/%0.4i\n',t,nt)
%         for p=1+2*xwidth*twidth:np-2*xwidth*twidth+2*nxplus
%             opts = optimset('Display','off');
%             flux(t,p) = fminbnd(@fastenergy4,-xwidth,xwidth,opts,yff,t,p,twidth,xwidth);
%         end
%     end
%     %fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')

    hackfact = 5;

    fu = zeros(nt, np); fv = zeros(nt,np);
    [yx yt] = gradient(yff); 
    yx = yx*hackfact;
    ytytf = imfilter(yt.*yt,h,'replicate');
    ytyxf = imfilter(yt.*yx,h,'replicate');
    yxyxf = imfilter(yx.*yx,h,'replicate');
    for k=1:nt*np
        A = [ytytf(k) ytyxf(k); ytyxf(k) yxyxf(k)];
        [u s v] = svd(A);
        fu(k) = u(1,1); fv(k) = u(2,1);
        fv(k) = fv(k) / hackfact;
    end
    
    edges(i).fu = fu;
    edges(i).fv = fv;
    
    if nargout>1
        x = edges(i).points2;
        %flux = repmat(-fu./fv,[1 1 2]) .* repmat(shiftdim(gradient(x),-1),[nt 1 1]);
        xind = round(x); xind = xind(:,2) + nj*(xind(:,1)-1);
        flux = -fu./fv;
        FLUX(:,xind) = flux;
    end
    
end
warning on MATLAB:divideByZero

if nargout>1
    FLUX = reshape(FLUX,nt,nj,ni);
end
    