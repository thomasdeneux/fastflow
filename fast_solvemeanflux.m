function [E, bins] = fast_solvemeanflux(e,stim,rest,sigmat,sigmax,SEflag,immean)
% function [E bins] = fast_solvemeanflux(e,stim,rest,sigmat,sigmax[,SEflag,immean])

nst = length(stim);
nre = length(rest);
[nt np] = size(e(1).data);

% flux moyen calculé de toute façon sur le global !
fprintf('[global] ')
E(1).ytyt = fn_means(e(stim).ytyt);
E(1).ytyx = fn_means(e(stim).ytyx);
E(1).yxyx = fn_means(e(stim).yxyx);
E(2).ytyt = fn_means(e(rest).ytyt);
E(2).ytyx = fn_means(e(rest).ytyx);
E(2).yxyx = fn_means(e(rest).yxyx);
E = fast_solve2D(E,sigmat,sigmax,'A2f');
E(3).flux = (E(1).flux - E(2).flux) / (mean(E(2).flux(:) + E(1).flux(:))/2);

% compute standard error ?
if nargin<7, immean = zeros(nt,length(e)); end
if nargin<6 || ~SEflag
    % volume
    E(1).volume = repmat(mean(immean(:,stim),2),1,np)+fn_means(e(stim).data);
    E(2).volume = repmat(mean(immean(:,rest),2),1,np)+fn_means(e(rest).data);
    E(3).volume = E(1).volume ./ E(2).volume;
    return
end

nbins = floor(min(nst,nre)/8); % nbins = 8 pour Exp1

% estimation du flux par bins
for k=1:nbins
    stk = stim(1+floor((k-1)*nst/nbins):floor(k*nst/nbins));
    rek = rest(1+floor((k-1)*nre/nbins):floor(k*nre/nbins));
    bins(k,1) = struct('ytyt',fn_means(e(stk).ytyt),'ytyx',fn_means(e(stk).ytyx),'yxyx',fn_means(e(stk).yxyx), ...
        'volume',repmat(mean(immean(:,stk),2),1,np)+fn_means(e(stk).data));
    bins(k,2) = struct('ytyt',fn_means(e(rek).ytyt),'ytyx',fn_means(e(rek).ytyx),'yxyx',fn_means(e(rek).yxyx), ...
        'volume',repmat(mean(immean(:,rek),2),1,np)+fn_means(e(rek).data));
end
fprintf('[bins] ')
bins = fast_solve2D(bins,sigmat,sigmax,'A2f');
for k=1:nbins
    bins(k,3).flux = (bins(k,1).flux - bins(k,2).flux) / (mean(bins(k,1).flux(:)+bins(k,2).flux(:))/2);
    bins(k,3).volume = bins(k,1).volume ./ bins(k,2).volume;
end

% erreurs
flux1 = reshape(cat(3,bins(:,1).flux),nt*np,nbins);
%E(1).flux = reshape(mean(flux1,2),nt,np);
E(1).fluxSE = reshape(sqrt(var(flux1')')/sqrt(nbins),nt,np);
flux0 = reshape(cat(3,bins(:,2).flux),nt*np,nbins);
% E(2).flux = reshape(mean(flux0,2),nt,np);
E(2).fluxSE = reshape(sqrt(var(flux0')')/sqrt(nbins),nt,np);
flux = reshape(cat(3,bins(:,3).flux),nt*np,nbins);
% E(3).flux = reshape(mean(flux,2),nt,np);
E(3).fluxSE = reshape(sqrt(var(flux')')/sqrt(nbins),nt,np);

volume1 = reshape(cat(3,bins(:,1).volume),nt*np,nbins);
E(1).volume = reshape(mean(volume1,2),nt,np);
E(1).volumeSE = reshape(sqrt(var(volume1')')/sqrt(nbins),nt,np);
volume0 = reshape(cat(3,bins(:,2).volume),nt*np,nbins);
E(2).volume = reshape(mean(volume0,2),nt,np);
E(2).volumeSE = reshape(sqrt(var(volume0')')/sqrt(nbins),nt,np);
volume = reshape(cat(3,bins(:,3).volume),nt*np,nbins);
E(3).volume = reshape(mean(volume,2),nt,np);
E(3).volumeSE = reshape(sqrt(var(volume')')/sqrt(nbins),nt,np);

