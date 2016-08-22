function Y = fast_anisosmooth(Yname,CSU,CSV,mask)
% function Y = fast_anisosmooth(Yname,CSU,CSV[,mask])
%----
% spatial smoothing of blood volume time courses data
% orthogonally to flow constraint
% Y must be space x space x time

% passage par r�f�rence �trange -> save memory !!
if ischar(Yname)
    Y = evalin('caller',Yname);
    evalin('caller',['clear ' Yname])
else
    Y = Yname;
end

% other inputs
if ~exist('mask','var')
    mask = [];
end
if ischar(mask), error('prototype has changed'), end

% filter design
filterhw = 3;
filterwidth = 2*filterhw+1;
h = zeros(filterwidth);
%h(filterhw+1,:)=fspecial('gaussian',[1 filterwidth],1);
h(filterhw+1,filterhw:filterhw+2)=1;
H = zeros(filterwidth,filterwidth,100);
[u v] = meshgrid(-filterhw:filterhw,-filterhw:filterhw);
u=u(:);v=v(:);
for i=1:101
    theta = pi*(i-1)/100;
    A = [cos(theta) -sin(theta) ; sin(theta) cos(theta)];
    uv2 = A*([u'; v']);
    hh = reshape(interp2(h,uv2(1,:)+filterhw+1,uv2(2,:)+filterhw+1),filterwidth,filterwidth);
    hh(isnan(hh))=0; H(:,:,i) = hh / sum(hh(:));
end

% Go !
[nj ni nt] = size(Y);
if ~isempty(mask), 
    Y = reshape(Y,nj*ni,nt);
    Y(find(~mask),:)=0;
    Y = reshape(Y,nj,ni,nt);
end
fprintf('lissage:           ')
for i=2+filterhw:ni-filterhw
    fprintf('\b\b\b\b\b\b\b\b\b\b%4d/%4d\n',i,ni)
    for j=1+filterhw:nj-filterhw
        if isempty(mask) || Y(j,i,1)~=0
            theta = round(100*atan(-CSV(j,i)/CSU(j,i))/pi+51);
            % H(:,:,theta) est un filtre gaussien orthogonal au vecteur
            % (CSU(j,i),CSV(j,i)), on filtre en faisant attention �
            % n'utiliser que les pixels dans le masque 

            % 3 following lines are replaced by MEX file call
            cross = repmat(H(:,:,theta),[1 1 nt]) .* Y(j-filterhw:j+filterhw,i-filterhw:i+filterhw,:);
            cross = reshape(cross,filterwidth^2,nt);
            Y(j,i-filterhw-1,:) = shiftdim(sum(cross),-1);
            %fast_anisosmooth_sub(Y,H(:,:,theta),j,i,filterhw)
            
            if ~isempty(mask)
                ref = H(:,:,theta) .* mask(j-filterhw:j+filterhw,i-filterhw:i+filterhw);
                ref = ref(:);
                Y(j,i-filterhw-1,:) = Y(j,i-filterhw-1,:) / sum(ref); 
            end
        end
    end
    pause(0)
end

fprintf('bascule:           ')
for i=ni-filterhw:-1:2+filterhw
    fprintf('\b\b\b\b\b\b\b\b\b\b%4d/%4d\n',i,ni)
    Y(:,i,:) = Y(:,i-filterhw-1,:);
end