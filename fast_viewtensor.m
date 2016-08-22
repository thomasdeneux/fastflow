function fast_viewtensor(edge,sub)
% function fast_viewtensor(edge,sigma,sub)

if nargin<2, sub=[10 10]; end
imagesc(edge.dataf)
hold on, fn_tensordisplay(edge,'sub',sub,'r'), hold off
