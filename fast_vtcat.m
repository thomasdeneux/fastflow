function e = fast_vt(vessels,trials,mode)
% function e = fast_vt(vessels,trials,'r'|'s'|'rs')

spath = [fn_cd('mars') '/save/catflow'];
i=0;

fn_progress('loading vessel:',length(vessels))
for w = vessels(:)'
    i=i+1;
    fn_progress(i)
    names = {};
    
    k = 1;
    if findstr(mode,'r')        
        for t=trials(:)'
            names{end+1} = ['v' num2str(w) '_t' num2str(t) '_r'];
        end
        load([spath '/edge' num2str(w,'%.2i')],names{:})
        j=0;
        
        np = length(eval(['v' num2str(w) '_t' num2str(trials(1)) '_r.points2']));
        for t=trials(:)'
            j=j+1;
            e(i,j,k).name = ['vessel' num2str(w) '_trial' num2str(t) '_r'];
            e(i,j,k).np=np;
            e(i,j,k).points = eval(['v' num2str(w) '_t' num2str(t) '_r.points']);
            e(i,j,k).points2 = eval(['v' num2str(w) '_t' num2str(t) '_r.points2']);
            e(i,j,k).data = eval(['v' num2str(w) '_t' num2str(t) '_r.data']);
            e(i,j,k).flux = eval(['v' num2str(w) '_t' num2str(t) '_r.flux']);
        end
        clear v*
        k = k+1;
    end
        
    if findstr(mode,'s')      
        for t=trials(:)'
            names{end+1} = ['v' num2str(w) '_t' num2str(t) '_s'];
        end
        load([spath '/edge' num2str(w,'%.2i')],names{:})
        j=0;
        
        np = length(eval(['v' num2str(w) '_t' num2str(trials(1)) '_s.points2']));
        for t=trials(:)'
            j=j+1;
            e(i,j,k).name = ['vessel' num2str(w) '_trial' num2str(t) '_s'];
            e(i,j,k).np=np;
            e(i,j,k).points = eval(['v' num2str(w) '_t' num2str(t) '_s.points']);
            e(i,j,k).points2 = eval(['v' num2str(w) '_t' num2str(t) '_s.points2']);
            e(i,j,k).data = eval(['v' num2str(w) '_t' num2str(t) '_s.data']);
            e(i,j,k).flux = eval(['v' num2str(w) '_t' num2str(t) '_s.flux']);
        end
        clear v*
    end
        
    end
end

e = fast_solve2D(e,30,5,'A');

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
