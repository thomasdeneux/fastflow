function fast_displayang(varargin)
% function fast_displayang([taxis,[img,]]ang)
% function fast_displayang([taxis,]res)
%---
% displays estimation result E obtained with fast_solvemeanflux

[nt np] = size(E(3).flux);
subt = 30:60:nt;

vol = 2-mean(E(3).volume,2);
volf = filty(vol,60,'lm');
volSE = filty(-mean(E(3).volumeSE,2),60,'lm');
flux = mean(E(3).flux,2);
fluxf = filty(flux,60,'lm');
fluxSE = filty(mean(E(3).fluxSE,2),60,'lm');

figure(1)
hold off, plot(taxis,vol,'color',[0 .5 0]), hold on
plot(taxis,volf,'linewidth',2,'color',[0 .5 0])
h = errorbar(taxis(subt),volf(subt),volSE(subt),'k');
delete(h(2))
axis tight, xlabel('time (s)'), ylabel('% signal change')
title('volume response to stimulation')

figure(2)
hold off, plot(taxis,flux), hold on
plot(taxis,fluxf,'linewidth',2)
h = errorbar(taxis(subt),fluxf(subt),fluxSE(subt),'k');
delete(h(2))
axis tight, xlabel('time (s)'), ylabel('% signal change')
title('flow response to stimulation')
