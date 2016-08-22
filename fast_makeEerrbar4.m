function [E, tot] = fast_makeEerrbar4(j,stim,rest,p,savepath)
% function E = fast_makeEerrbar(j,stim,rest,p,savepath)
%---
% Input:
% - j               vessel number
% - stim, rest      indices of trials (stimulated and rest)
% - p               parameters for structure tensor smoothing
% - savepath        directory wher the files 'marsedgeIVOnew_v??_t??.mat' are
% 
% Output:
% - E               1x3 structure with fields ytyt, ytyx, yxyx,
%                   flux,fluxtot,errflux

nst = length(stim);
nre = length(rest);
nbins = floor(min(nst,nre)/8); % nbins = 8 pour Exp1

bins = cell(3,nbins); % (flow stim, flow rest, flow stim/rest) for each bin
tot = struct('ytyt',{0 0},'ytyx',{0 0},'yxyx',{0 0});

% estimation du flux par bins
for k=1:nbins
    disp(['bin ' num2str(k) '/' num2str(nbins)])
    stk = stim(1+floor((k-1)*nst/nbins):floor(k*nst/nbins));
    rek = rest(1+floor((k-1)*nre/nbins):floor(k*nre/nbins));
    
    % average flow for stimulated trials inside bin k
    ebin = struct('ytyt',0,'ytyx',0,'yxyx',0);
    fn_progress('reading stim',length(stk))
    for ii=1:length(stk)
        fn_progress(ii)
        i = stk(ii);
        % load
        vtnameload=['v' num2str(j) '_t' num2str(i)];
%         if i < 205
%             load([savepath '\edges1\marsedgeIVOnew' '_v' num2str(j,'%.2i') '_t' num2str(i,'%.3i')],vtnameload)
%         else
%             load([savepath '\edges2\marsedgeIVOnew' '_v' num2str(j,'%.2i') '_t' num2str(i,'%.3i')],vtnameload)
%         end
        load([savepath '\marsedgeIVOnew' '_v' num2str(j,'%.2i') '_t' num2str(i,'%.3i')],vtnameload)
        
        e = eval(vtnameload); clear(vtnameload)
        
        % sum 
        if i < 205      
            ebin.ytyt = ebin.ytyt + e.ytyt(91:1200,:)/length(stk);        
            ebin.ytyx = ebin.ytyx + e.ytyx(91:1200,:)/length(stk);        
            ebin.yxyx = ebin.yxyx + e.yxyx(91:1200,:)/length(stk);        
            tot(1).ytyt = tot(1).ytyt + e.ytyt(91:1200,:)/length(stim);        
            tot(1).ytyx = tot(1).ytyx + e.ytyx(91:1200,:)/length(stim);        
            tot(1).yxyx = tot(1).yxyx + e.yxyx(91:1200,:)/length(stim);
        else
            ebin.ytyt = ebin.ytyt + e.ytyt(1:1110,:)/length(stk);        
            ebin.ytyx = ebin.ytyx + e.ytyx(1:1110,:)/length(stk);        
            ebin.yxyx = ebin.yxyx + e.yxyx(1:1110,:)/length(stk);        
            tot(1).ytyt = tot(1).ytyt + e.ytyt(1:1110,:)/length(stim);        
            tot(1).ytyx = tot(1).ytyx + e.ytyx(1:1110,:)/length(stim);        
            tot(1).yxyx = tot(1).yxyx + e.yxyx(1:1110,:)/length(stim);
        end
    end
     
    ebin = fast_dotensornodata(ebin,p,'A2f');
    bins{1,k} = ebin.flux;
    
    % average flow for rest trials inside bin k
    ebin = struct('ytyt',0,'ytyx',0,'yxyx',0);
    fn_progress('reading rest',length(rek))
    for ii =1:length(rek)
        fn_progress(ii)
        i = rek(ii);
        % load
        vtnameload=['v' num2str(j) '_t' num2str(i)];
%         if i < 205
%             load([savepath '\edges1\marsedgeIVOnew' '_v' num2str(j,'%.2i') '_t' num2str(i,'%.3i')],vtnameload)
%         else
%             load([savepath '\edges2\marsedgeIVOnew' '_v' num2str(j,'%.2i') '_t' num2str(i,'%.3i')],vtnameload)
%         end
        load([savepath '\marsedgeIVOnew' '_v' num2str(j,'%.2i') '_t' num2str(i,'%.3i')],vtnameload)
        e = eval(vtnameload); clear(vtnameload)
        % sum 
        if i < 205
            ebin.ytyt = ebin.ytyt + e.ytyt(91:1200,:)/length(rek);
            ebin.ytyx = ebin.ytyx + e.ytyx(91:1200,:)/length(rek);       
            ebin.yxyx = ebin.yxyx + e.yxyx(91:1200,:)/length(rek);        
            tot(2).ytyt = tot(2).ytyt + e.ytyt(91:1200,:)/length(rest);        
            tot(2).ytyx = tot(2).ytyx + e.ytyx(91:1200,:)/length(rest);        
            tot(2).yxyx = tot(2).yxyx + e.yxyx(91:1200,:)/length(rest);
        else
            ebin.ytyt = ebin.ytyt + e.ytyt(1:1110,:)/length(rek);
            ebin.ytyx = ebin.ytyx + e.ytyx(1:1110,:)/length(rek);       
            ebin.yxyx = ebin.yxyx + e.yxyx(1:1110,:)/length(rek);        
            tot(2).ytyt = tot(2).ytyt + e.ytyt(1:1110,:)/length(rest);        
            tot(2).ytyx = tot(2).ytyx + e.ytyx(1:1110,:)/length(rest);        
            tot(2).yxyx = tot(2).yxyx + e.yxyx(1:1110,:)/length(rest);
        end
           
    end
    ebin = fast_dotensornodata(ebin,p,'A2f');
    bins{2,k} = ebin.flux;
        
    % average sensory-evoked flow inside bin k
    bins{3,k} = (bins{1,k}-bins{2,k}) ; %/ mean(bin(1,k)+bin(2,k));
end
    
% average and standard error over bins
disp('averaging over bins')
[nt np] = size(bins{1,1});
flux1 = reshape(cat(3,bins{1,:}),nt*np,nbins);
E(1).flux = reshape(mean(flux1,2),nt,np);
E(1).fluxSE = reshape(sqrt(var(flux1')')/sqrt(nbins),nt,np);
flux0 = reshape(cat(3,bins{2,:}),nt*np,nbins);
E(2).flux = reshape(mean(flux0,2),nt,np);
E(2).fluxSE = reshape(sqrt(var(flux0')')/sqrt(nbins),nt,np);
flux = reshape(cat(3,bins{3,:}),nt*np,nbins);
E(3).flux = reshape(mean(flux,2),nt,np);
E(3).fluxSE = reshape(sqrt(var(flux')')/sqrt(nbins),nt,np);

% % si le flux estimé n'est pas bon, essayer plutot un moyennage sur tous les
% % trials
 tot = fast_dotensornodata(tot,p,'A2f'); 
 E(1).fluxtot = tot(1).flux;
 E(2).fluxtot = tot(2).flux;
 E(3).fluxtot = tot(1).flux-tot(2).flux;
