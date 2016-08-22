function e = fast_vt(vessels,trials,flag)
% function e = fast_vt(vessels,trials,'older'|'old'|'new')

if nargin<3, flag = 'old'; end
switch flag
    case 'older'
        edgefilename = [fn_cd('mars') '/save/8nov_fastflow/edge'];
    case 'old'
        edgefilename = [fn_cd('mars') '/save/5dec_fastflow/edge'];
    case 'new'
        edgefilename = ['E:\Users\Thomas\0503_Marseille\Exp1\janvedge'];
end
i=0;
fn_progress('loading vessel:',length(vessels))
for w = vessels(:)'
    i=i+1;
    fn_progress(i)
        
    if strcmp(flag,'aug')
        basename = ['E:\Users\Thomas\0503_Marseille\'];
        j=0;
        for t=trials(:)'
            j=j+1;
            e(i,j).name = ['vessel' num2str(w) '_trial' num2str(t)];
            load([basename 'ivo_edges\marsedgeIVOnew_v' num2str(w,'%.2i') '_t' num2str(t,'%.3i')])
            e(i,j).data = eval(['v' num2str(w) '_t' num2str(t) '.data']);
            load([basename 'ivo_edges_tensor\marsedgeIVOnew_v' num2str(w,'%.2i') '_t' num2str(t,'%.3i')])
            e(i,j).ytyt = eval(['v' num2str(w) '_t' num2str(t) '.ytyt']);
            e(i,j).ytyx = eval(['v' num2str(w) '_t' num2str(t) '.ytyx']);
            e(i,j).yxyx = eval(['v' num2str(w) '_t' num2str(t) '.yxyx']);
            clear v*
        end
    else
        names = {};
        for t=trials(:)'
            names{end+1} = ['v' num2str(w) '_t' num2str(t)];
        end
        load([edgefilename num2str(w,'%.2i')],names{:})
        j=0;
        np = length(eval(['v' num2str(w) '_t' num2str(trials(1)) '.points2']));
        for t=trials(:)'
            j=j+1;
            e(i,j).name = ['vessel' num2str(w) '_trial' num2str(t)];
            e(i,j).np=np;
            e(i,j).points = eval(['v' num2str(w) '_t' num2str(t) '.points']);
            e(i,j).points2 = eval(['v' num2str(w) '_t' num2str(t) '.points2']);
            e(i,j).data = eval(['v' num2str(w) '_t' num2str(t) '.data']);
            try, e(i,j).flux = eval(['v' num2str(w) '_t' num2str(t) '.flux']); end
        end
        clear v*
    end
end

%e = fast_solve2D(e,30,5,'A');

