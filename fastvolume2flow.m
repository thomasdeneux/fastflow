function FLUX = fastvolume2flow(Yf,csu,csv,alpha,beta,gamma)
% function FLUX = fastvolume2flow(Y,csu,csv,alpha,beta,gamma)
%---
% DESUETE, use fast_2Dsolve

% Estimation of blood speed from blood volume movie
% Input:
% - Y               3-D array (time x space x space)
% - csu, csv        2-D array defining the directions constraint on flow estimation
% - alpha, beta     energy coefficient for fit to data and flow smoothness,

[nt ny nx] = size(Yf);

% % visualization
% figure(1), clf reset
% set(1,'doublebuffer','on','Name','filtered movie')
% colormap gray
% mM = [min(Yf(:)) max(Yf(:))];
% for i=1:nt, imagesc(squeeze(Yf(i,:,:)),mM), pause(0.03), end
% pause

% % derivatives
% % (spécifique fastenergy)
% [Yy Yt Yx] = gradient(Yf);
% v0x = repmat(shiftdim(csu,-1),nt,1);
% v0y = repmat(shiftdim(csv,-1),nt,1);
% Yv0 = Yx.*v0x + Yy.*v0y;

% optical flow : contrainte v = k v0 (le long de la vasculature)
% dI/dx . v + dI/dt = 0

% % résolution directe: k = - dI/Dt / (dI/dx . v0)
% FLUX = - Yt ./ Yv0;dd

% descente de gradient avec continuité spatiale
% % (spécifique fastenergy 1/2)
% flux0 = zeros(ny,nx);
% FLUX = zeros(nt,ny,nx);
flux0 = zeros(ny,nx,2);
FLUX = zeros(nt,ny,nx,2);

% % options de minimisation
% % (spécifique fastenergy)
% jacobpattern = sparse(2*nx*ny,nx*ny);
% for i=1:nx
%     for j=1:ny
%         jacobpattern(j+ny*(i-1),j+ny*(i-1))=1;
%         for ii=max(1,i-1):min(nx,i+1)
%             for jj=max(1,j-1):min(ny,j+1)
%                  jacobpattern((nx*ny)+jj+ny*(ii-1),j+ny*(i-1))=1;
%             end
%         end
%     end
% end
% opts = optimset('JacobPattern',jacobpattern,'Display','off');

% options de minimisation
% (spécifique fastenergy2)
opts = optimset('Jacobian','on','Display','off');

disp('flow estimation          ')
for i=1:nt-1
    fprintf('\b\b\b\b\b\b\b\b\b\b%0.4d/%0.4d\n',i,nt)
% % (spécifique fastenergy)
%     yt = squeeze(Yt(i,:,:));
%     yv0 = squeeze(Yv0(i,:,:));
%     flux = lsqnonlin(@fastenergy,flux0,[],[],opts,yv0,yt,csu,csv,alpha,beta);
    flux = lsqnonlin(@fastenergy3,flux0,[],[],opts, ...
        squeeze(Yf(i,:,:)),squeeze(Yf(i+1,:,:)),csu,csv,alpha,beta,gamma);
    FLUX(i,:,:,:) = shiftdim(flux,-1);
end

save tmp20septFLUX FLUX

% % visualization
% figure(1), clf reset
% set(1,'Name','Structural constraints for flow')
% colormap gray
% fn_4Dview('quiver',1,csu,csv,squeeze(Yf(1,:,:)))
% fn_4Dview('4D',2,3,Yf);
% figure(2), colormap gray
% set(2,'Name','Volume map')
% set(3,'Name','Volume time course')
% fn_4Dview('activeplot',4,FLUX)
% figure(4), colormap gray
% set(4,'Name','Flow time course')
% fn_4Dview('quiver',5,FLUX,csu,csv,squeeze(Yf(1,:,:)))
% set(5,'Name','Flow map')
% pause
% 

