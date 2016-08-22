function mask = fast_morphomask(mask,nconx,nclose)
% function mask = fast_morphomask(mask[,nconx,nclose])
%---
% remove small connex components and performs closure to fill holes

if nargin==0, help fast_morphomask, end

if nargin<2, nconx=50; end
if nargin<3, nclose=5; end

if nconx
    [lab ncomp] = bwlabel(mask);
    compt = zeros(1,666);
    for k=lab(:)', if k, compt(k) = compt(k)+1; end, end
    badlab = find(compt<nconx);
    mask(ismember(lab,badlab))=0;
end

for i=1:nclose
    mask = bwmorph(mask,'dilate');
end
for i=1:nclose
    mask = bwmorph(mask,'erode');
end
