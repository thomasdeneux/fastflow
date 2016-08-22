function Y = mik_registration(Y)
% function Y = mik_registration(Y)

disp('save rawfloat movie')
fast_saverawfloat(Y,'tmp.rawfloat')

disp('sf_registration tmp.rawfloat')
!sf_registration tmp.rawfloat

disp('read rawfloat movie')
Y = fast_loadrawfloat('tmp-filtered.rawfloat');
%delete tmp.rawfloat tmp.match tmp-registered.rawfloat tmp-filtered.rawfloat