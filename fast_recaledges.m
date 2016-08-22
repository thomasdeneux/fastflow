function edges = fast_recaledges(edges,iso,nj,ni)
% function edges = fast_recaledges(edges,iso,nj,ni)
% function edges = fast_recaledges(edges,M)

if nargin<4
    M = iso;
else 
    M = fast_register(iso,nj,ni);
end

for i=1:length(edges)
    edges(i).points = [edges(i).points ones(size(edges(i).points,1),1)]*M';
    edges(i).points2 = [edges(i).points2 ones(size(edges(i).points2,1),1)]*M';
end
