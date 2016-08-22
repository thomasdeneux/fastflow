function avrg = TOM_avrgresults
% function avrg = TOM_avrgresults

anduillette=pwd; %#ok<NASGU>
reslist=myls('res_TOM_*.mat');
% reslist = reslist(2:2:end,:);

y=load(reslist(1,:));y=y.res;
nsc = length(y.anglesc_syl);

anglesc=y.anglesc;mult_anglesc=(y.anglesc~=0);
anglesc_orig=y.anglesc_orig;mult_anglesc_orig=(y.anglesc_orig~=0);
speedsc=y.speedsc;mult_speedsc=(y.speedsc~=0);
angletg=y.angletg;mult_angletg=(y.angletg~=0);
speedtg=y.speedtg;mult_speedtg=(y.speedtg~=0);
for k = 1:nsc
    anglesc_syl{k}=y.anglesc_syl{k};mult_anglesc_syl{k}=(y.anglesc_syl{k}~=0); %#ok<*AGROW>
    speedsc_syl{k}=y.speedsc_syl{k};mult_speedsc_syl{k}=(y.speedsc_syl{k}~=0);
    angletg_syl{k}=y.angletg_syl{k};mult_angletg_syl{k}=(y.angletg_syl{k}~=0);
    speedtg_syl{k}=y.speedtg_syl{k};mult_speedtg_syl{k}=(y.speedtg_syl{k}~=0);
end
fbest=(y.bestang*180/pi~=1); % don't remember: bestang=pi where estimation failed?
bestmult=fbest;
bestang=zeros(size(y.bestang)); bestspeed=zeros(size(y.bestang));
bestang(fbest)=y.bestang(fbest);
bestspeed(fbest)=tan(y.bestang(fbest));

for i = 2:size(reslist,1)
    y=load(reslist(i,:));y=y.res;
    anglesc=(anglesc+y.anglesc);mult_anglesc=mult_anglesc+(y.anglesc~=0);
    anglesc_orig=(anglesc_orig+y.anglesc_orig);mult_anglesc_orig=mult_anglesc_orig+(y.anglesc_orig~=0);
    speedsc=(speedsc+y.speedsc);mult_speedsc=mult_speedsc+(y.speedsc~=0);
    angletg=(angletg+y.angletg);mult_angletg=mult_angletg+(y.angletg~=0);
    speedtg=(speedtg+y.speedtg);mult_speedtg=mult_speedtg+(y.speedtg~=0);
    for k = 1:4
        anglesc_syl{k}=(anglesc_syl{k}+y.anglesc_syl{k});mult_anglesc_syl{k}=mult_anglesc_syl{k}+(y.anglesc_syl{k}~=0);
        speedsc_syl{k}=(speedsc_syl{k}+y.speedsc_syl{k});mult_speedsc_syl{k}=mult_speedsc_syl{k}+(y.speedsc_syl{k}~=0);
        angletg_syl{k}=(angletg_syl{k}+y.angletg_syl{k});mult_angletg_syl{k}=mult_angletg_syl{k}+(y.angletg_syl{k}~=0);
        speedtg_syl{k}=(speedtg_syl{k}+y.speedtg_syl{k});mult_speedtg_syl{k}=mult_speedtg_syl{k}+(y.speedtg_syl{k}~=0);
    end
    fbest=(y.bestang*180/pi~=1);
    bestmult=bestmult+fbest;
    bestang(fbest)=bestang(fbest)+y.bestang(fbest);
    bestspeed(fbest)=bestspeed(fbest)+tan(y.bestang(fbest));
end

avrg.anglesc=anglesc./(mult_anglesc+isinf(ones(size(mult_anglesc))./mult_anglesc));
avrg.anglesc_orig=anglesc_orig./(mult_anglesc_orig+isinf(ones(size(mult_anglesc_orig))./mult_anglesc_orig));
avrg.speedsc=speedsc./(mult_speedsc+isinf(ones(size(mult_speedsc))./mult_speedsc));
avrg.angletg=angletg./(mult_angletg+isinf(ones(size(mult_angletg))./mult_angletg));
avrg.speedtg=speedtg./(mult_speedtg+isinf(ones(size(mult_speedtg))./mult_speedtg));
for k = 1:4
    avrg.anglesc_syl{k}=anglesc_syl{k}./(mult_anglesc_syl{k}+1e-6);
    avrg.speedsc_syl{k}=speedsc_syl{k}./(mult_speedsc_syl{k}+1e-6);
    avrg.angletg_syl{k}=angletg_syl{k}./(mult_angletg_syl{k}+1e-6);
    avrg.speedtg_syl{k}=speedtg_syl{k}./(mult_speedtg_syl{k}+1e-6);
end
avrg.bestang=bestang./(bestmult + (bestmult==0)) * 180/pi;
avrg.bestspeed=bestspeed./(bestmult + (bestmult==0));

% figure
% y=plot(compute_speed(avrg.anglesc,'angle')*0.8);set(y,'color','k','Linewidth',2,'Linestyle','-')
% hold
% y=plot(compute_speed(avrg.anglesc,'angle')*0.8);set(y,'color','k','Linewidth',2,'Linestyle','-')
% y=plot(compute_speed(avrg.bestang,'angle')*0.8);set(y,'color','b','Linewidth',2,'Linestyle','-')
% grid
% title(['angles- ',anduillette(49:end)])
% legend('angsc','angsc orig',1);
