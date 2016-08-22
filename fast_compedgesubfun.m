function fast_displayedgesubfun(j)

FLUX = evalin('base','FLUX');
FLUX2 = evalin('base','FLUX2');
% 
% cd C:\DATA\newfastflow\anl
% 
% %j = evalin('caller','j');
% 
% load diversTOTAL17mar06
% Edgepath='C:\DATA\newfastflow\edgesE12\'
% %timeaxisnew=[1:2399]*0.0025-0.5;
timeaxisnew=[1:2399-180]*0.0025-0.05;
% cd(Edgepath)
% fileslist=myls('E*');
nedges=size(FLUX,1);

% load(fileslist(j+1,:))

figure(2), clf
subplot(1,2,1)
yy=plot(timeaxisnew,FLUX(j,1).flux*.8,'r');set(yy,'Linewidth',2);
hold
yy=plot(timeaxisnew,FLUX(j,2).flux*.8,'b');set(yy,'Linewidth',2);
yy=plot(timeaxisnew,FLUX(j,1).fluxtot*.8,'r--');set(yy,'Linewidth',2);
yy=plot(timeaxisnew,FLUX(j,2).fluxtot*.8,'b--');set(yy,'Linewidth',2);
xlabel('time after stim. onset(s)')
ylabel('RBC speed (pix/frame)')
title(sprintf(['vessel n. ' num2str(j) '; r=stim; b=rest; --=totflux']))
axis tight

subplot(1,2,2)
fluxbin=FLUX(j,3).flux*.8;
dummy = mean(FLUX(j,1).flux + FLUX(j,2).flux)*.8/2;
yy=plot(timeaxisnew,fluxbin/dummy,'m-');set(yy,'Linewidth',2);
hold
fluxtot=FLUX(j,3).fluxtot*.8;
dummytot = mean(FLUX(j,1).fluxtot + FLUX(j,2).fluxtot)*.8/2;
yy=plot(timeaxisnew,fluxtot/dummytot,'r');set(yy,'Linewidth',2);
SE=FLUX(j,3).fluxSE*.8;
yy=plot(timeaxisnew,fluxbin/dummy-SE/abs(dummy),'g-');set(yy,'Linewidth',2);
yy=plot(timeaxisnew,fluxbin/dummy+SE/abs(dummy),'g-');set(yy,'Linewidth',2);
yy=plot(timeaxisnew,fluxtot/dummytot-SE/abs(dummytot),'c-');set(yy,'Linewidth',2);
yy=plot(timeaxisnew,fluxtot/dummytot+SE/abs(dummytot),'c-');set(yy,'Linewidth',2);
lowsig=100;
yy=plot(timeaxisnew,filty(fluxbin/dummy,lowsig,'lm'),'m-');set(yy,'Linewidth',2);
yy=plot(timeaxisnew,filty(fluxtot/dummytot,lowsig,'lm'),'r');set(yy,'Linewidth',2);
yy=plot(timeaxisnew,filty(fluxbin/dummy-SE/abs(dummy),lowsig,'lm'),'g-');set(yy,'Linewidth',2);
yy=plot(timeaxisnew,filty(fluxbin/dummy+SE/abs(dummy),lowsig,'lm'),'g-');set(yy,'Linewidth',2);
yy=plot(timeaxisnew,filty(fluxtot/dummytot-SE/abs(dummytot),lowsig,'lm'),'c-');set(yy,'Linewidth',2);
yy=plot(timeaxisnew,filty(fluxtot/dummytot+SE/abs(dummytot),lowsig,'lm'),'c-');set(yy,'Linewidth',2);

xlabel('time after stim. onset(s)')
ylabel('RBC speed response (pix/frame)')
title(sprintf(['vessel n. ' num2str(j) '; STIM-REST']))
axis tight

timeaxisnew=[1:1110]*0.005-0.05;

figure(3), clf
subplot(1,2,1)
yy=plot(timeaxisnew,FLUX2(j,1).flux*.8,'r');set(yy,'Linewidth',2);
hold
yy=plot(timeaxisnew,FLUX2(j,2).flux*.8,'b');set(yy,'Linewidth',2);
yy=plot(timeaxisnew,FLUX2(j,1).fluxtot*.8,'r--');set(yy,'Linewidth',2);
yy=plot(timeaxisnew,FLUX2(j,2).fluxtot*.8,'b--');set(yy,'Linewidth',2);
xlabel('time after stim. onset(s)')
ylabel('RBC speed (pix/frame)')
title(sprintf(['vessel n. ' num2str(j) '; r=stim; b=rest; --=totflux']))
axis tight

subplot(1,2,2)
fluxbin=FLUX2(j,3).flux*.8;
dummy = mean(FLUX2(j,1).flux + FLUX2(j,2).flux)*.8/2;
yy=plot(timeaxisnew,fluxbin/dummy,'m-');set(yy,'Linewidth',2);
hold
fluxtot=FLUX2(j,3).fluxtot*.8;
dummytot = mean(FLUX2(j,1).fluxtot + FLUX2(j,2).fluxtot)*.8/2;
yy=plot(timeaxisnew,fluxtot/dummytot,'r');set(yy,'Linewidth',2);
SE=FLUX2(j,3).fluxSE*.8;
yy=plot(timeaxisnew,fluxbin/dummy-SE/abs(dummy),'g-');set(yy,'Linewidth',2);
yy=plot(timeaxisnew,fluxbin/dummy+SE/abs(dummy),'g-');set(yy,'Linewidth',2);
yy=plot(timeaxisnew,fluxtot/dummytot-SE/abs(dummytot),'c-');set(yy,'Linewidth',2);
yy=plot(timeaxisnew,fluxtot/dummytot+SE/abs(dummytot),'c-');set(yy,'Linewidth',2);
lowsig=100;
yy=plot(timeaxisnew,filty(fluxbin/dummy,lowsig,'lm'),'m-');set(yy,'Linewidth',2);
yy=plot(timeaxisnew,filty(fluxtot/dummytot,lowsig,'lm'),'r');set(yy,'Linewidth',2);
yy=plot(timeaxisnew,filty(fluxbin/dummy-SE/abs(dummy),lowsig,'lm'),'g-');set(yy,'Linewidth',2);
yy=plot(timeaxisnew,filty(fluxbin/dummy+SE/abs(dummy),lowsig,'lm'),'g-');set(yy,'Linewidth',2);
yy=plot(timeaxisnew,filty(fluxtot/dummytot-SE/abs(dummytot),lowsig,'lm'),'c-');set(yy,'Linewidth',2);
yy=plot(timeaxisnew,filty(fluxtot/dummytot+SE/abs(dummytot),lowsig,'lm'),'c-');set(yy,'Linewidth',2);

xlabel('time after stim. onset(s)')
ylabel('RBC speed response (pix/frame)')
title(sprintf(['vessel n. ' num2str(j) '; STIM-REST']))
axis tight
