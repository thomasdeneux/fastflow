function FLUX = fast_flowest(Y,CSU,CSV,mask) 
% function FLUX = fast_flowest(Y[,CSU,CSV[,mask]]) 
%---
% flow estimation from fast volume acquisition
% ABANDONNED - use fast_2Dsolve instead

% Input:
% - Y           blood volume data (time x space x space)
% - CSU, CSV    flow direction constraint (space x space)
%               alternatively, use a single (space x space x 2) array
%               estimated if not given
% Output:
% - flux        blood flow vector field (time x space x space x 2)


% essai ultime (?) - ABANDONNé !!!!
% on n'extrait plus de lignes 1D
% en chaque point, on crée un 1D local en utilisant l'information de
% contrainte
% ensuite on estime le flux en utilisant la matrice de structure (?)

% parameters
streamstep = 1.5;
streamlength = 5;
ns = 2*streamlength-1;
thighcut = 15;
tlowcut = 1;
xlowcut = 1;
h = fspecial('gaussian',[8 8],3); % flow smoothing

% Input
[nt ny nx] = size(Y);
Y = permute(Y,[2 3 1]);
if max(max(Y(:,:,1)))>3 || min(min(Y(:,:,1)))<.3
    % flow direction constraint
    if nargin<2
        [CSU CSV] = fast_structure(Y(:,:,1));
    elseif nargin<3
        CSV = CSU(:,:,2);
        CSU = CSU(:,:,1);
    end
    % normalisation    
    Y = Y ./ repmat(mean(Y),nt,1);
elseif nargin<3
    CSV = CSU(:,:,2);
    CSU = CSU(:,:,1);
end

% problem: CSU and CSV may not be continous
% -> define 2 possible new continous vector fields
[CSUH CSUV] = deal(CSU);
[CSVH CSVV] = deal(CSV);
f = find(CSU<0); CSUH(f)=-CSU(f); CSVH(f)=-CSV(f);
f = find(CSV<0); CSUV(f)=-CSU(f); CSVV(f)=-CSV(f);

FLUX = zeros(nt,ny,nx,1);
% ATTENTION gérer mieux les sorties de range pour éviter le 20 ci-dessous
disp('estimation          ')
for y=20:ny-20
    fprintf('\b\b\b\b\b\b\b\b\b\b%4i/%4i\n',y,ny)
    for x=20:nx-20
        
        if nargin>=4 && ~mask(y,x), continue, end
        
        % data along vessel
        if abs(CSU(y,x))>abs(CSV(y,x))
            CSUa = CSUH; CSVa = CSVH;
        else
            CSUa = CSUV; CSVa = CSVV;
        end
        strm1 = stream2(CSUa,CSVa,x,y,[streamstep streamlength]);
        strm2 = stream2(-CSUa,-CSVa,x,y,[streamstep streamlength]);
        % ATTENTION gérer mieux la longueur des segments !
        stream = [strm2{1}(end:-1:2,:) ; strm1{1}];
%         yy = zeros(nt,ns);
%         for t=1:nt
%             yy(t,:) = interp2(Y(:,:,t),stream(:,1),stream(:,2))';
%         end
        yy = tmpinterp2(Y,stream)';
        
        % smoothing
        yf = filtx(filty(filty2(yy(:,:),thighcut,'hm'),tlowcut,'lm'),xlowcut,'lm');

        % flow estimation
        [fu fv] = deal(zeros(nt,1));
        [ys yt] = gradient(yf);
        ytytf = imfilter(yt.*yt,h,'replicate');
        ytysf = imfilter(yt.*ys,h,'replicate');
        ysysf = imfilter(ys.*ys,h,'replicate');
        for t=1:nt
            A = [ytytf(t,1+streamlength) ytysf(t,1+streamlength); ...
                    ytysf(t,1+streamlength) ysysf(t,1+streamlength)];
            [u s v] = svd(A);
            fu(t) = u(1,1); fv(t) = u(2,1);
            FLUX(t,y,x) = -u(1,1)/u(2,1);
        end
        
    end
end

