function [F, DF] = fastenergy2(FLUX,IA,IB,CSU,CSV,alpha,beta,gamma)
% function [F DF] = fastenergy2(FLUX,IA,IB,CSU,CSV,alpha,beta,gamma)
%---
% Energie et dérivée basée sur:
% - fit to data: cross-correlation entre IA(x) et IB(x+v(x))
% - smoothness:  divergence du flux

[nj ni] = size(IA);
[I J] = meshgrid(1:ni,1:nj);
VI = FLUX(:,:,1);
VJ = FLUX(:,:,2);

% IB(x+v(x))
IB2 = interp2(IB,I+VI,J+VJ);

% correlation
Fflow = IB2-IA;

% smoothness: divergence
Fsmooth = divergence(VI,VJ);

% direction constraint: vectoriel product
Fdir = VI.*CSV - VJ.*CSU;


if nargout>1
    delta = 1e-3;
    IB2b = interp2(IB,I+VI+delta,J+VJ);
    IB2X = (IB2b-IB2)/delta;
    IB2b = interp2(IB,I+VI,J+VJ+delta);
    IB2Y = (IB2b-IB2)/delta;
   
    % correlation
    DFflow = spdiags([IB2X(:) IB2Y(:)],[0 ni*nj],ni*nj,2*ni*nj);
    
    % smoothness
    DFsmooth = sparse([],[],[],nj*ni,2*nj*ni,4*ni*nj);
    for j=1:nj
        for i=1:ni
            k = j+nj*(i-1);
            k2 = k+nj*ni;
            if j==1
                DFsmooth(k,k2) = -1;
                DFsmooth(k,k2+1) = 1;
            elseif j==nj
                DFsmooth(k,k2-1) = -1;
                DFsmooth(k,k2) = 1;
            else
                DFsmooth(k,k2-1) = -1/2;
                DFsmooth(k,k2+1) = 1/2;
            end    
            if i==1
                DFsmooth(k,k) = DFsmooth(k,k)-1;
                DFsmooth(k,k+nj) = 1;
            elseif i==ni
                DFsmooth(k,k-nj) = -1;
                DFsmooth(k,k) = DFsmooth(k,k)+1;
            else
                DFsmooth(k,k-nj) = -1/2;
                DFsmooth(k,k+nj) = 1/2;
            end    
        end
    end
    
    % direction
    DFdir = spdiags([CSV(:) -CSU(:)],[0 ni*nj],ni*nj,2*ni*nj);
    
    
end

% corrections
%f = isnan(Fflow(:));
Fflow([1:4 end-3:end],:)=0;
Fflow(:,[1:4 end-3:end])=0;
F = [alpha*Fflow(:) ; beta*Fsmooth(:) ; gamma*Fdir(:)];
if nargout>1
    [indi indj] = meshgrid([1:4 ni-3:ni],1:nj); ind = indj + nj*(indi-1);
    DFflow(ind,:)=0;
    [indi indj] = meshgrid(1:ni,[1:4 nj-3:nj]); ind = indj + nj*(indi-1);
    DFflow(ind,:)=0;
    DF = [alpha*DFflow ; beta*DFsmooth ; gamma*DFdir];
end

% % à supprimer
% F = reshape(full(F),nj,ni,3);
% if nargout>1
%     DF = reshape(full(DF),nj,ni,3,nj,ni,2);
% end


