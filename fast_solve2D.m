function edges = fast_solve2D(edges,sigmat,sigmax,flag)
% function edges = fast_solve2D(edges[,sigmat,sigmax[,'A'|'f'|'A2f']])
%---
% Blood flow estimation from fast blood volume measure
% 
% Input:
% - edges   structure with field 'data' containing time courses in edges
% - sigmat  temporal smoothing
% - sigmax  spatial smoothing
% - flag    string
%
% Output:
% - edges   fields 'fu', 'fv' (estimated "directions" in 'data')
%           are added (or 'A' is flag is set to 'A')


% flow smoothing (more exactly: represents the extent of neighbourghood
% when making structure matrice, see below)
%  time x space
if nargin<3
    sigmat = 30;
    sigmax = 5;
end
if nargin<4
    flag = 'f';
end
h = fspecial('gaussian',[ceil(2*sigmat) 1],sigmat)*fspecial('gaussian',[1 ceil(2*sigmax)],sigmax);

disp('attention: spatial filtering has changed')
disp('estimating flow           ')
nedges = prod(size(edges));

warning off MATLAB:divideByZero
for i=1:nedges
    fprintf('\b\b\b\b\b\b\b\b\b\b%4i/%4i\n',i,nedges)
    switch flag
        case {'A','f'}
            yf = edges(i).data;    
            [nt np] = size(yf);
            
            yff = imfilter(filty2(yf(:,:),35,'hm'),fspecial('gaussian',[5 1],1)*fspecial('gaussian',[1 5],1));
            %yff = imfilter(filty2(yf(:,:),35,'hm'),fspecial('gaussian',[5 1],1));
            %yff = filt2(filty2(yf(:,:),35,'hm'),.5,'lm');
            %yff = filty2(yf(:,:),35,'hm');
            
            [yx yt] = gradient(yff); 
            yx = yx;
            ytyt = yt.*yt;
            ytyx = yt.*yx;
            yxyx = yx.*yx;
        case 'A2f'
            ytyt = edges(i).ytyt;
            ytyx = edges(i).ytyx;
            yxyx = edges(i).yxyx;
            [nt np] = size(ytyt);
    end            
    switch flag
        case {'f','A2f'}
            hackfact = 1;
            ytytf = imfilter(ytyt,h,'replicate');
            ytyxf = imfilter(ytyx,h,'replicate')*hackfact;
            yxyxf = imfilter(yxyx,h,'replicate')*(hackfact^2);
            fu = zeros(nt, np); fv = zeros(nt,np); %poids = zeros(nt,np);
            for k=1:nt*np
                A = [ytytf(k) ytyxf(k); ytyxf(k) yxyxf(k)];
                [u s v] = svd(A);
                poids = s(1,1)/s(2,2);
                fu(k) = u(1,1); fv(k) = u(2,1);
                fv(k) = fv(k) / hackfact;
            end
            norme = sqrt(fu.^2 + fv.^2);
            if strcmp(flag,'A2f'), edges(i).fu = fu./norme; end
            if strcmp(flag,'A2f'), edges(i).fv = fv./norme; end
            edges(i).flux = -fu./fv;
            
%             % � supprimer + tard
%             edges(i).fu = fu./norme; 
%             edges(i).fv = fv./norme; 
%             edges(i).datafiltered = yff;
            
        case 'A'
            edges(i).ytyt = ytyt;
            edges(i).ytyx = ytyx;
            edges(i).yxyx = yxyx;
    end    
end
warning on MATLAB:divideByZero

% opts = optimset('Display','off');
% for i=1:nedges
%     fprintf('\b\b\b\b\b\b\b\b\b\b%4i/%4i\n',i,nedges)
%     yf = edges(i).data;
%     yff = imfilter(filty2(yf(:,:),15,'hm'),fspecial('gaussian',[5 1],1)*fspecial('gaussian',[1 5],1));
%     
%     np = size(yff,2);
% 
%     flux = zeros(nt, 1);
%     for t=1:nt-1
%         M = Inf; ind = 0;
%         y1 = yff(t,:); y2 = yff(t+1,:);
%         for dec = -10:10
%             E = lsqenergy(dec,y1,y2);
%             if E<M, M=E; ind=dec; end
%         end
%         flux(t) = ind;
%     end
%     edges(i).flux = flux;
%         
%     if nargout>1
%         x = edges(i).points2;
%         %flux = repmat(-fu./fv,[1 1 2]) .* repmat(shiftdim(gradient(x),-1),[nt 1 1]);
%         xind = round(x); xind = xind(:,2) + nj*(xind(:,1)-1);
%         FLUX(:,xind) = repmat(flux,1,np);
%     end
% end


%-------------------------------------------------------------------------
function E=lsqenergy(dec,y1,y2)

np = length(y1);
if dec<0, extrapval=y1(1); else extrapval=y1(end); end
y1dec = interp1(1:np,y1,(1:np)+dec,'cubic',extrapval);
E = norm(y2(:)-y1dec(:));