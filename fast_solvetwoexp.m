function [E, bins] = fast_solvetwoexp(e1,e2,stim1,stim2,rest1,rest2,sigmat,sigmax,immean)
% function [E, bins] = fast_solvetwoexp(e1,e2,stim1,stim2,rest1,rest2,sigmat,sigmax,immean)
%---
% Input:
% - e1              n1 struct
% - e2              n2 struct
% - stim1, rest1    indices in [1 n1]
% - stim2, rest2    indices in [1 n2]
% - sigmat, sigmax  scalars
% - immean          nt x (n1+n2) array

cut1 = 91:1200;
cut2 = 1:1110;

n1 = length(e1);
for i=1:n1
    e1(i).ytyt = e1(i).ytyt(cut1,:);
    e1(i).ytyx = e1(i).ytyx(cut1,:);
    e1(i).yxyx = e1(i).yxyx(cut1,:);
    e1(i).data = e1(i).data(cut1,:);
end
immean1 = immean(cut1,1:n1);

n2 = length(e2);
for i=2:n2
    e2(i).ytyt = e2(i).ytyt(cut2,:);
    e2(i).ytyx = e2(i).ytyx(cut2,:);
    e2(i).yxyx = e2(i).yxyx(cut2,:);
    e2(i).data = e2(i).data(cut2,:);
end
immean2 = immean(cut2,n1+1:n1+n2);

stim = [stim1(:) ; stim2(:)+n1];
rest = [rest1(:) ; rest2(:)+n1];
immean = [immean1 immean2];

e = [e1 e2];
E = fast_solvemeanflux(e,stim,rest,sigmat,sigmax,true,immean);
