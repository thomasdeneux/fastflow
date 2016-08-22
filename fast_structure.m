function [mask, skel, nodes] = fast_structure(CS,mask,skel)
% function [maskbin skel nodes] = fast_structure(CS[,mask[,skel]])
% function [maskbin CSU CSV] = fast_structure(CS,'csu')
%--------------------------------
% performs work on a static vasculature image:
% - vessels orientation
% - vessels segmentation
%
% Input:
% - CS          2D image (first dimension = y)
%
% Output:
% - CSU,CSV     vector fiel = constraint for blood flow direction

if nargin<1, help fast_structure, return, end

% sizes
[nj ni] = size(CS);

%
% VESSELS ORIENTATION
%

if nargin>=2 && ischar(mask)
    flag = lower(mask);
    clear mask
else
    flag = '';
end

if ~exist('mask','var')
    disp('vessels orientation')
    [I J] = meshgrid(1:ni,1:nj);
    step = 1;
    [tmpi tmpj] = meshgrid(1:step:ni,1:step:nj); sub = tmpj + nj*(tmpi-1);
    
    % smoothing
    h = fspecial('gaussian',[5 5],.6);
    CS = imfilter(CS,h,'replicate');
    CS = CS / mean(CS(:));
    
    % derivatives
    [CSI CSJ] = gradient(CS);
    [CSII CSIJ] = gradient(CSI); %#ok<NASGU>
    [CSIJ CSJJ] = gradient(CSJ);
    
    % vessels directions
    CSU = CS; CSV = CS;
    for i=1:ni*nj
        A = [CSII(i) CSIJ(i); CSIJ(i) CSJJ(i)];
        [u s v] = svd(A);
        CSU(i) = u(1,2); CSV(i) = u(2,2);
    end
    
    % vasculature maskbin
    mask = filt2(CSII+CSJJ,1);
end

maskbin = (mask>0);

if nargout<2 || strcmp(flag,'csu')
    skel = CSU;
    edges = CSV;
    return
end

%
% VESSELS SKELETTON -> GRAPH
%

if ~exist('skel','var')
    
    % remove small connex components
    disp('performing binary operations')
    [lab ncomp] = bwlabel(maskbin);
    count = hist(lab(:),0:ncomp);
    [jnk ord] = sort(count(2:end));
    maskbin = ismember(lab,ord(end-10:end));
    
    if nargout<2, return, end
    
    skel0 = bwmorph(maskbin,'skel',Inf);
    
    % eliminating small branches
    skel = elague(skel0,5);
    
end

if nargout<3, return, end

disp('find extremities')
lut = makelut(@lut_design,3,'+3');
cross = applylut(skel,lut);
lut = makelut(@lut_design,3,'1');
ends  = applylut(skel,lut);
nodes = cross | ends;
nodes = bwmorph(nodes,'shrink',Inf);
[ii jj] = find(nodes);
nodes = [ii jj];


% % find cross-sections = graph nodes
% disp('labelling edges')
% lut = makelut(@lut_design,3,'+3');
% cross = applylut(skel,lut);
% 
% % labelling graph edges
% pieces = xor(skel,cross); % cross-sections removed
% [lab nedges] = bwlabel(pieces); %#ok<NASGU>
% 
% % remove short edges
% count = hist(lab(:),0:max(lab(:)));
% count(1)=[];
% f = find(count>10);
% pieces = ismember(lab,f); % short components removed
% % re-label
% newval = zeros(1,nedges+1);
% nedges = length(f);
% newval(f+1) = 1:nedges;
% lab = newval(lab+1);
% 
% % sorting according to labels
% lut = makelut(@lut_design,3,'-1');
% ends = applylut(pieces,lut);
% [yends xends] = find(ends);
% edgends = zeros(nedges,2);
% edgends(lab(ends),:) = [xends yends]; % extr�mit�s tri�es par labels
% % qd plusieurs ext/edge, la 2�me recouvre la 1�re
% 
% % collecting successive points on edges
% disp('collecting edges')
% edges = struct('points',cell(nedges,1),'points2',[],'ext',[0 0]);
% dirs = [1 1 1 0 -1 -1 -1 0
%     -1 0 1 1 1 0 -1 -1]';
% for k=1:nedges
%     npoints = count(k);
%     xs = zeros(npoints,2);
%     xs(1,:) = edgends(k,:);
%     inddir0 = 0; % direction dans laquelle chercher en  dernier, 
%     % �vitera de retourner en arri�re par la suite
%     for j=2:npoints
%         x = xs(j-1,:);
%         for inddir = 1+mod(inddir0+(0:7),8)
%             xd = x+dirs(inddir,:);
%             if any(xd<1 | xd>[ni nj]), continue, end
%             if pieces(xd(2),xd(1)), break, end
%         end
%         xs(j,:) = xd;
%         inddir0 = mod(inddir+4,8);
%     end
%     edges(k).points = xs;
% end
% 
% % connecting nodes
% disp('collecting nodes')
% [ycross xcross] = find(cross);
% nnodes = length(xcross);
% nodes = struct('coords',mat2cell([xcross ycross],ones(nnodes,1),2), ...
%     'edg',[],'edgdir',logical([]),'vois',[]);
% for i=1:nnodes
%     xy = nodes(i).coords;
%     voisedges = voisins(pieces,xy); n1 = size(voisedges,1);
%     voiscross = voisins(cross,xy); n2 = size(voiscross,1);
%     edg = zeros(1,n1);
%     edgdir = zeros(1,n1);
%     vois = zeros(1,n2);
%     for j=1:n1
%         xyv = voisedges(j,:);
%         edg(j) = lab(xyv(2),xyv(1));
%         edgdir(j) = all(edges(edg(j)).points(1,:)==xyv);
%         if size(edges(edg(j)).points,1)==1 && edges(edg(j)).ext(1)~=0
%             % cas particulier: edge n'a qu'un point qui repr�sente les 2
%             % extr�mit�s !
%             edges(edg(j)).ext(2)=i;
%         else
%             edges(edg(j)).ext(2-edgdir(j)) = i;
%         end
%     end
%     for j=1:n2
%         xyv = voiscross(j,:);
%         vois(j) = find(xcross==xyv(1) & ycross==xyv(2));
%     end
%     nodes(i).edg = edg;
%     nodes(i).edgdir = edgdir;
%     nodes(i).vois = vois;
% end
% 
% % making nice edges
% disp('resample edges')
% for k=1:nedges
%     %x = [nodes(edges(k).ext(1)).coords ; edges(k).points ; nodes(edges(k).ext(2)).coords];
%     x = edges(k).points;
%     x = nicesnake(x,3.5);
%     x = nicesnake(x,.7);
%     edges(k).points2 = x;
% end


%-------------------------
function x2 = nicesnake(x,ds)

np = size(x,1);
if np==1, x2=x; return, end
L = zeros(np-1,1);
for i=2:np, L(i) = L(i-1)+norm(x(i,:)-x(i-1,:)); end
if ~isempty(L), x2 = interp1(L,x,0:ds:L(end),'spline'); end
if norm(x(end,:)-x2(end,:))<ds/3
    x2(end,:)=x(end,:);
else
    x2 = [x2 ; x(end,:)];
end

%-------------------------
function skel = elague(skel,n)
% remove small branches from skeletton skel whos length are <= n

keep = false(size(skel));

shrink = bwmorph(skel,'shrink',1);
for i=1:n
    rem = skel - shrink;
    skel = shrink;
    shrink = bwmorph(skel,'shrink',1);
    % we want to keep the parts which belong to 'long branches', i.e. which
    % are disconnected from skel once skel has been 'shrinked' one step
    % further
    lab = bwlabel(shrink | rem);
    keep(~ismember(lab,[0; lab(shrink)])) = 1;
end
% only the 'real long branches', i.e. parts of 'keep' which connect to
% 'skel' are to be kept
lab = bwlabel(skel | keep);
skel(ismember(lab,lab(skel))) = 1;

%-------------------------
function vois = voisins(I,x,y)
% find neighbourghs of x in a logical array
% convention image: dimension order = y, x

if nargin>=3, x=[x y]; else x=x(:)'; end

vois = [];
for d = [1 1 1 0 -1 -1 -1 0 ; -1 0 1 1 1 0 -1 -1]
    xd = x+d';
    try
        if I(xd(2),xd(1)), vois = [vois ; xd]; end
    end
end

