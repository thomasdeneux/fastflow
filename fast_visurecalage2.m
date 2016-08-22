% script fast_visurecalage
% - doit etre exécuté dans le répertoire où sont les fichiers divers.mat et recalage??.mat
% - calcule une variable 'mask' qui est le masque de la surface du cerveau
%   qui reste visible tout au long de l'expérience 
% - affiche les paramètres de recalage de tous les trials dans une fenetre
%   . on peut zoomer...
%   . quand on clique sur la courbe, le film du trial correspondant est
%     chargé et on visualise avec le recalage
%     ensuite une fenetre demande si on veut supprimer ce trial ou non
%   . une variable 'badtrials' contient les noms des trials qu'on veut supprimer

% load 
load divers
iso = zeros(nt,3,nexp2);
for i=1:nexp2
    load(sprintf('recalagebis%.2i',ord2(i)))
    iso(:,:,i)=ISO;
end
clear ISO

% mask
fprintf('computing mask           ')
% subax = [1 ni 1 nj];
% for i=1:nexp
%     fprintf('\b\b\b\b\b\b\b\b%3i/%3i\n',i,nexp)
%     for j=1:nt
%         M = fast_register(iso(j,:,i),nj,ni);
%         M(3,:) = [0 0 1]; M = M^-1;
%         carre = M*[1 1 ni ni ; 1 nj 1 nj ; 1 1 1 1];
%         subax(1) = max([subax(1) ceil(carre(1,1:2))]);
%         subax(2) = min([subax(2) floor(carre(1,3:4))]);
%         subax(3) = max([subax(3) ceil(carre(2,[1 3]))]);
%         subax(4) = min([subax(4) floor(carre(2,[2 4]))]);
%     end
% end
iso = permute(iso,[1 3 2]);
iso = reshape(iso,nt*nexp2,3);
iso(1,:) = [0 0 0];
subax = [1-floor(min(iso(:,2))) ni-ceil(max(iso(:,2))) 1-floor(min(iso(:,3))) nj-ceil(max(iso(:,3)))];
mask = false(nj,ni);
mask(subax(3):subax(4),subax(1):subax(2)) = true;

% display
figure(2), colormap gray, set(2,'doublebuffer','on')
figure(1)
fn_imvalue
taxis = (0:nt*nexp2-1)/nt +.5;
hl=plot(taxis,[iso(:,1)*500 iso(:,2:3)]);
legend('rotation','translation x','translation y')
axis tight, ax = axis;
for i=0:nexp
    line([i i]+.5,ax(3:4),'color','black')
end
clear ax
if ~exist('Y','var') || ndims(Y)~=3 || ~all(size(Y)==[nj ni nt])
    Y = zeros(nj,ni,nt);
end
bad = [];
set(hl,'HitTest','on','ButtonDownFcn', [...
    'p=get(gca,''currentpoint'');' ...
    'k=round(p(1));' ...
    'i=ord2(round(p(1)));' ...
    'answer = questdlg(''remove trial ?'',num2str(i),''Yes'',''No'',''Play movie'',''No'');' ...
    'if strcmp(answer,''Play movie'')' ...
    '   Y=fast_loaddata_ref(''Y'',files2(i,:),1:nt);' ...
    '   Y=fast_recalage(''Y'',iso((1:nt)+nt*(k-1),:));' ...
    '   figure(2);' ...
    '   answer = ''Replay'';' ...
    '   while strcmp(answer,''Replay''),' ...
    '      for j=1:nt, imagesc(Y(:,:,j)), pause(0), end;' ...
    '      answer = questdlg(''remove trial ?'','''',''Yes'',''No'',''Replay'',''No'');' ...
    '   end,' ...
    'end,' ...
    'if strcmp(answer,''Yes''), bad=union(bad,i); end,' ...
])




