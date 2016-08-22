function e = fast_doall(e,fname)
% function e = fast_doall(e,fname)
%---
% Input:
% - e       structure with fields points (not needed actually), points2
% - fname   name of data file
%
% Output:
% - e       fields are added:
%           data (data interpolation)
%           

swd = pwd;
cd \\chicoree\deneux\WORK\IRMf\0503_Marseille\data\fastflow
disp('find indices')
CS = fast_loaddata(fname,1);
[nj ni] = size(CS);
[e indices] = fast_edgesmouse(e,nj,ni);
disp('load data')
Y = fast_loaddata(fname,800,indices,'normalize');
e = fast_interpalongedges(e,Y,indices,nj,ni);
e = fast_solve2D(e,30,5,'A');
e = fast_solve2D(e,30,5,'A2f');
cd(swd)