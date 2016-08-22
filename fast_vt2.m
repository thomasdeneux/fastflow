function e = fast_vt2(vessels,trials)
% function e = fast_vt2(vessels[,trials])

spath = [fn_cd('mars') '/save/5dec_fastflow'];
i=0;
for w = vessels(:)'
    i=i+1;
    names = {};
    for t=trials(:)'
        names{end+1} = ['v' num2str(w) '_t' num2str(t)];
    end
    load([spath '/edgebis' num2str(w,'%.2i')],names{:})
    j=0;
    %np = Inf;
    %for t=trials(:)'
    %    np = min(length(eval(['v' num2str(w) '_t' num2str(t) '.points2'])),np);
    %end
    np = length(eval(['v' num2str(w) '_t' num2str(trials(1)) '.points2']));
    for t=trials(:)'
        j=j+1;
        e(i,j).name = ['vessel' num2str(w) '_trial' num2str(t)];
        e(i,j).np=np;
        e(i,j).points = eval(['v' num2str(w) '_t' num2str(t) '.points']);
        e(i,j).points2 = eval(['v' num2str(w) '_t' num2str(t) '.points2']);
        e(i,j).data = eval(['v' num2str(w) '_t' num2str(t) '.data']);
        e(i,j).flux = eval(['v' num2str(w) '_t' num2str(t) '.flux']);
    end
    e = fast_solve2D(e,30,5,'A');
    clear v*
end

%     np = Inf;
%     for t=trials(:)'
%         np = min(np,size(eval(['v' num2str(w) '_t' num2str(t) '.points2']),1));
%     end
%     j=0;
%     fnames = fieldnames(eval(['v' num2str(w) '_t' num2str(t)]));
%     for t=trials(:)'
%         j=j+1;
%         e(i,j).name = ['vessel' num2str(w) '_trial' num2str(t)];
%         e(i,j).points = eval(['v' num2str(w) '_t' num2str(t) '.points']);
%         e(i,j).points2 = eval(['v' num2str(w) '_t' num2str(t) '.points2(1:np,:)']);
%         e(i,j).data = eval(['v' num2str(w) '_t' num2str(t) '.data(:,1:np)']);
%         e(i,j).ytyt = eval(['v' num2str(w) '_t' num2str(t) '.ytyt(:,1:np)']);
%         e(i,j).ytyx = eval(['v' num2str(w) '_t' num2str(t) '.ytyx(:,1:np)']);
%         e(i,j).yxyx = eval(['v' num2str(w) '_t' num2str(t) '.yxyx(:,1:np)']);
%         e(i,j).flux = eval(['v' num2str(w) '_t' num2str(t) '.flux(:,1:np)']);
%     end
