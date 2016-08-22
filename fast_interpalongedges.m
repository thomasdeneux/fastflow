function [edges] = fast_interpalongedges(edges,Y,ISO)
% function [edges] = fast_interpalongedges(edges,Y[,ISO])
%---
% Blood flow estimation from fast blood volume measure
% 
% Input:
% - Y       blood volume (space x space x time)
% - edges   structure with field 'points2' containing segmentation of vessels
%           it is estimated if not specified           
% - ISO     rotation + translation parameters of Y registering (will be
%           applied to edges points)
%
% Output:
% - edges   'data' field is added (interpolated data)

if nargin==0, help fast_interpalongedges, return, end

nedges = length(edges);

[nj ni nt] = size(Y);
mask = (Y(:,:,1)~=0);

warningedges=[];

% data
edges(1).data=[];
%disp('interpolating edges          ')
for i=1:nedges
    %fprintf('\b\b\b\b\b\b\b\b\b\b%4i/%4i\n',i,nedges)
    if any(edges(i).points2(:,1)<3 | edges(i).points2(:,1)>ni-2 | ...
            edges(i).points2(:,2)<3 | edges(i).points2(:,2)>nj-2)
        warningedges = [warningedges i];
        edges(i).points2(:,1) = min(max(edges(i).points2(:,1),3),ni-2);
        edges(i).points2(:,2) = min(max(edges(i).points2(:,2),3),nj-2);
    end
    x = edges(i).points2;
    nx = size(x,1);
    if nargin<3
        data = zeros(nt,nx);
        for k=1:nx
            qi = floor(x(k,1)); ri = x(k,1)-qi;
            qj = floor(x(k,2)); rj = x(k,2)-qj;
            f = find(mask(qj:qj+1,qi:qi+1));
            Y0 = Y(qj:qj+1,qi:qi+1,:);
            Y0 = permute(Y0,[3 1 2]);
            r = [(1-rj) rj]'*[(1-ri) ri]; 
            if sum(r(f))==0
                r(:)=nan;
            else
                r = r/sum(r(f));
            end
            for z=f'
                data(:,k) = data(:,k) + r(z)*Y0(:,z);
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
    else % on ne g�re pas le masque dans ce cas
        % x2 est 2 x nx x nt
        x2 = fast_recalage(x,ISO,nj,ni);
        % x2qi, ..., indices sont nx x nt
        x2q = floor(x2); x2qi = shiftdim(x2q(1,:,:),1); x2qj = shiftdim(x2q(2,:,:),1);
        x2r = x2-x2q; x2ri = shiftdim(x2r(1,:,:),1); x2rj = shiftdim(x2r(2,:,:),1);
        clear x2q x2r
        indices = x2qj + nj*(x2qi-1) + nj*ni*repmat([0:nt-1],nx,1); 
        % data est nt x nx
        try
            data = (1-x2ri) .* ((1-x2rj).*Y(indices) + x2rj.*Y(indices+1)) ...
                + x2ri .* ((1-x2rj).*Y(indices+nj) + x2rj.*Y(indices+nj+1));
        catch 
            data = zeros(nx,nt);
            disp(['probl�me: le vaisseau ' num2str(i) ' sort du champ de l''image !'])
        end
        data = data';
    end
    edges(i).data = data;
end

if ~isempty(warningedges), warning(['edges ' num2str(warningedges,'%i ') 'were too close to border']), end