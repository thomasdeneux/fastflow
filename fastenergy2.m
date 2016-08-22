function [F, DF] = fastenergy2(FLUX,IA,IB,CSU,CSV,alpha,beta)
% function [F DF] = fastenergy2(FLUX,IA,IB,CSU,CSV,alpha,beta)
%---
% Energie et dérivée basée sur:
% - fit to data: cross-correlation entre IA(x) et IB(x+v(x))
% - smoothness:  divergence du flux

[nj ni] = size(IA);
[I J] = meshgrid(1:ni,1:nj);
VI = FLUX.*CSU;
VJ = FLUX.*CSV;

% IB(x+v(x))
IB2 = interp2(IB,I+VI,J+VJ);

% % cross-correlation 1
% % IAxIB2(x) = A([IA(voisx) ; IB(voisx)]) 
% %           = IA(voisx) - IB(voisx) - mean(IA(voisx) - IB(voisx))
% nvois = 3; % doit etre impair 
% nvois2 = nvois^2;
% kvois = (nvois-1)/2;
% A = [eye(nvois2) -eye(nvois2)] ...
%     ; %- repmat([ones(1,nvois2) -ones(1,nvois2)]/(2*nvois2),nvois2,1);
% Fflow = zeros(nj,ni,nvois2);
% for j=nvois+2:nj-nvois-1
%     for i=nvois+2:ni-nvois-1
%         IAvoisx = IA(j-kvois:j+kvois,i-kvois:i+kvois);
%         IB2voisx = IB2(j-kvois:j+kvois,i-kvois:i+kvois);
%         Fflow(j,i,:) = shiftdim(A*[IAvoisx(:) ; IB2voisx(:)],-2);
%     end
% end

% % cross-correlation 2
% % IAxIB2(x) =  (dIA/dx - dIB2/dx).u(x)
% [IAX IAY] = gradient(IA);
% [IB2X IB2Y] = gradient(IB2);
% Fflow = (IAX-IB2X).*CSU + (IAY-IB2Y).*CSV;

% correlation tout court merde !
Fflow = IB2-IA;

% smoothness: divergence
Fsmooth = divergence(FLUX.*CSU,FLUX.*CSV);

if nargout>1
    delta = 1e-3;
    IB2b = interp2(IB,I+VI+delta*CSU,J+VJ+delta*CSV);
    DIB2 = (IB2b-IB2)/delta;

%     % cross-correlation 1
%     DFflow = sparse([],[],[],nvois2*nj*ni,nj*ni,nvois2*nvois2*nj*ni);
%     A2 = A(:,nvois2+1:end);
%     for j=nvois+2:nj-nvois-1
%         for i=nvois+2:ni-nvois-1
%             k = j+nj*(i-1);
%             % indices dans DFflow: correspondent à (j,i,pq,j+p,i+q)
%             % DFflow(j,i,:,j+jj,i+ii) = A2(:,:,jj,ii) * DIB2voisx(j+jj,i+ii,j+jj,i+ii)
%             for ii=-kvois:kvois
%                 for jj=-kvois:kvois
%                     DFflow(j+nj*((i-1)+ni*(0:nvois2-1)),(j+jj)+nj*(i+ii-1)) = ...
%                         A2(:,(jj+kvois+1)+nvois*(ii+kvois)) * DIB2(j+jj,i+ii);
%                 end
%             end
%         end
%     end

%     % cross-correlation 2
%     DFflow = sparse([],[],[],nj*ni,nj*ni,4*nj*ni);
%      for j=1:nj
%         for i=1:ni
%             k = j+nj*(i-1);
%             if j==1
%                 DFflow(k,k) = CSV(j,i) * DIB2(j,i);
%                 DFflow(k,k+1) = -CSV(j,i) * DIB2(j+1,i);
%             elseif j==nj
%                 DFflow(k,k-1) = CSV(j,i) * DIB2(j-1,i);
%                 DFflow(k,k) = -CSV(j,i) * DIB2(j,i);
%             else
%                 DFflow(k,k-1) = CSV(j,i)/2 * DIB2(j-1,i);
%                 DFflow(k,k+1) = -CSV(j,i)/2 * DIB2(j+1,i);
%             end    
%             if i==1
%                 DFflow(k,k) = DFflow(k,k) + CSU(j,i) * DIB2(j,i);
%                 DFflow(k,k+nj) = -CSU(j,i) * DIB2(j,i+1);
%             elseif i==ni
%                 DFflow(k,k-nj) = CSU(j,i) * DIB2(j,i-1);
%                 DFflow(k,k) = DFflow(k,k) - CSU(j,i) * DIB2(j,i);
%             else
%                 DFflow(k,k-nj) = CSU(j,i)/2 * DIB2(j,i-1);
%                 DFflow(k,k+nj) = -CSU(j,i)/2 * DIB2(j,i+1);
%             end    
%         end
%     end
   
    % correlation tout court merde !
    DFflow = sparse([],[],[],nj*ni,nj*ni,nj*ni);
    for k=1:nj*ni, DFflow(k,k)=DIB2(k); end
    
    DFsmooth = sparse([],[],[],nj*ni,nj*ni,4*ni*nj);
    for j=1:nj
        for i=1:ni
            k = j+nj*(i-1);
            if j==1
                DFsmooth(k,k) = -CSV(j,i);
                DFsmooth(k,k+1) = CSV(j+1,i);
            elseif j==nj
                DFsmooth(k,k-1) = -CSV(j-1,i);
                DFsmooth(k,k) = CSV(j,i);
            else
                DFsmooth(k,k-1) = -CSV(j-1,i)/2;
                DFsmooth(k,k+1) = CSV(j+1,i)/2;
            end    
            if i==1
                DFsmooth(k,k) = DFsmooth(k,k)-CSU(j,i);
                DFsmooth(k,k+nj) = CSU(j,i+1);
            elseif i==ni
                DFsmooth(k,k-nj) = -CSU(j,i-1);
                DFsmooth(k,k) = DFsmooth(k,k)+CSU(j,i);
            else
                DFsmooth(k,k-nj) = -CSU(j,i-1)/2;
                DFsmooth(k,k+nj) = CSU(j,i+1)/2;
            end    
        end
    end
    
    
end

% corrections
%f = isnan(Fflow(:));
Fflow([1:4 end-3:end],:)=0;
Fflow(:,[1:4 end-3:end])=0;
F = [alpha*Fflow(:) ; beta*Fsmooth(:)];
if nargout>1
    [indi indj] = meshgrid([1:4 ni-3:ni],1:nj); ind = indj + nj*(indi-1);
    DFflow(ind,:)=0;
    [indi indj] = meshgrid(1:ni,[1:4 nj-3:nj]); ind = indj + nj*(indi-1);
    DFflow(ind,:)=0;
    DF = [alpha*DFflow ; beta*DFsmooth];
end

% % à supprimer
% F = reshape(full(F),nj,ni,2);
% if nargout>1
%     DF = reshape(full(DF),nj,ni,2,nj,ni);
%     DF = permute(DF,[1 2 4 5 3]);
% end
% 

