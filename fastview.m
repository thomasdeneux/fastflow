function a = fastview(b)

% parameters structure
S = pointer(struct('step',''));

% menu figure
hf = figure(827);
m = fn_menu(hf,'v','FASTFLOW',50,20);
fn_menu(m,'add','string','NEW','callback',@(u,evnt)new(S))
fn_menu(m,'add','string','OPEN','callback',@(u,evnt)open(S))
fn_menu(m,'add','string','EDGES','callback',@(u,evnt)edges(S))
fn_menu(m,'add','string','COREGISTER','callback',@(u,evnt)coregister(S))
fn_menu(m,'add','string','ANALYSIS','callback',@(u,evnt)analysis(S))
    
% few initializations
set(0,'defaultfigurecolormap',gray(64))
fn_imvalue image

%---
function c=new(S)

S.x = struct('step','');
inputexperiment(S)
        
%---
function d=open(S)
            
fname = fn_getfile('*_batch.mat','Select existing batch data file');
if ~fname, return, end
load(fname)

% % 2 - data and saving directory
% 
% if strcmp(S.step,'')
%     S.datadir = fn_getfile('DIR',[],'Select data directory');
%     if strcmp(S.datadir,0), disp('batch interrupted'), return, end
%     S.datadir(end+1) = '/';
%     S.savedir = fn_getfile('DIR',[],'Select save directory');
%     if strcmp(S.savedir,0), disp('batch interrupted'), return, end
%     S.savedir(end+1) = '/';
%     
%     [dum base] = fileparts(S.datadir);
%     answer = inputdlg('Enter base nickname','',1,{base});
%     S.base = answer{1};
%     S.step = 'folders';
%     save([S.savedir S.base '_batch'],'S')
% end
% 
% % 3 - data files, conditions, sort according to experiments and conditions
% 
% if strcmp(S.step,'folders')
%     % data files
%     disp('file names')
%     a = dir(S.datadir);
%     f = false(1,length(a));
%     for k=1:length(a)
%         check = regexp(a(k).name,'^TC\d{2}\d{6}\d{6}_E\d{2}B\d{3}.BLK$');
%         f(k) = ~isempty(check);
%     end
%     
%     % select which experiments to take
%     files = strvcat(a(f).name);
%     exps = str2num(files(:,19:20))'; %#ok<ST2NM>
%     exp  = unique(exps);
%     answer = inputdlg('Select experiments','',1,{num2str(exp)});
%     exp = str2num(answer{1});
%     f = ismember(exps,exp);
%     files = files(f,:);
%     
%     % sort files
%     exps   = str2num(files(:,19:20)); %#ok<ST2NM>
%     blocks = str2num(files(:,22:24)); %#ok<ST2NM>
%     [dum ord] = sortrows([exps blocks]);
%     disp('Ordered file names')
%     disp(files(ord,:))
%     answer = questdlg('Approve files order?','','Yes','No','No');
%     if ~strcmp(answer,'Yes'), disp('batch interrupted'), return, end
%     S.files = files(ord,:);
%     S.nexp = size(S.files,1);
%     
%     % conditions
%     conds  = str2num(files(:,3:4));   %#ok<ST2NM>
%     cond = unique(conds)';
%     condstim = cond(2:end-1);
%     condrest = cond([1 end]);
%     answer = inputdlg({'stim conditions:','rest conditions:'}, ...
%         'Define conditions', ...
%         1, ...
%         {num2str(condstim),num2str(condrest)});
%     condstim = str2num(answer{1});
%     condrest = str2num(answer{2});
%     S.stim = find(ismember(conds,condstim));
%     S.rest = find(ismember(conds,condrest));
%     
%     S.step = 'files';
%     save([S.savedir S.base '_batch'],'S')
% end
% 
% % 4 - basic data
% 
% if strcmp(S.step,'files')
%     fname = [S.datadir S.files(1,:)];
%     hdr = fast_loaddata(fname,'header');
%     S.nx = hdr.xs;
%     S.ny = hdr.ys;
%     S.nt = hdr.nfrms;
%     S.CS = fast_loaddata(fname,1);
%     
%     [S.CSU S.CSV] = fast_structure(S.CS,'csu');
%     figure(1), imagesc(S.CS')
%     figure(2), fn_displayarrows(S.CS',{S.CSV' S.CSU'})
%     
%     % CHECK PROBLEM CSU IS LOGICAL !!!
%     
%     answer = questdlg('Approve reference frame?','','Yes','No','No');
%     if ~strcmp(answer,'Yes'), disp('batch interrupted'), return, end
%     
%     S.step = 'CS';
%     save([S.savedir S.base '_batch'],'S')
% end
% 
% % 5 - select vessels - to improve
% 
% if strcmp(S.step,'CS') || strcmp(S.step,'vessels definition')
%     if strcmp(S.step,'CS')
%         disp('process first trial')
%         Y = fast_loaddata([S.datadir S.files(1,:)]);
%         [iso Y] = fast_recalage(Y,S.CS);
%         S.tmp = struct('iso1',iso);
%         S.step = 'vessels definition';
%     elseif ~exist('Y','var')
%         disp('process first trial')
%         Y = fast_loaddata([S.datadir S.files(1,:)]);
%         Y = fast_recalage(Y,S.tmp.iso1);
%         X = fast_normalizeseq(Y,1,fspecial('gaussian',[1 11],3));
%     end
%     
%     disp('vessel selection:')
%     disp('press ''SAVE VESSEL'' to save the current vessel')
%     disp('press ''DONE'' when finished')
%     figure(1), clf
%     clip = [.996 1.004];
%     pos = [
%         .05 .05 .5  .45
%         .05 .53 .5  .45
%         .6 .05  .35 .9];
%     
%     
%     A = fn4D(axes('pos',pos(1,:)),Y,'2d'    ,'clip',[]);
%     C = fn4D(axes('pos',pos(2,:)),X,'2d'    ,'clip',clip);
%     D = fn4D(axes('pos',pos(3,:)),X,'snake' ,'clip',clip);
%     C.D.clip = clip; D.D.clip = clip;
%     
%     SI = D.SI;
%     
%     if ~exist('ke','var'), ke=0; end
%     uicontrol('units','normalized','position',[.91 .05 .07 .03],'string','SAVE VESSEL', ...
%         'callback','ke=ke+1; S.edges(ke) = struct(''snake'',SI.snake,''tube'',SI.tube,''A'',SI.A); save edg s ke');
%     u = uicontrol('units','normalized','position',[.91 .02 .07 .03],'string','DONE', ...
%         'callback',@(hu,evnt)delete(hu));
%     
% end
% % >> for k=1:8, S.edges(k).np = size(S.edges(k).snake,2); end
% % >> for k=1:8, A = S.edges(k).A; S.edges(k).a = sparse(A.p,A.i+S.nx*(A.j-1),A.z,S.edges(k).np,S.nx*S.ny); end
% 
% % 6 - coregister
% 
% if strcmp(S.step,'CS') || strcmp(S.step,'registering')
%     if strcmp(S.step,'CS')
%         S.ISO = zeros(S.nt,3,S.nexp);
%         kstart = 1;
%         iso = zeros(1,3);
%         S.step = 'registering';
%     else
%         kstart = find(~any(any(S.ISO)),1,'first');
%         iso = S.ISO(end,:,kstart-1);
%         answer = questdlg(['Continue interrupted registration? (trial ' num2str(kstart) ')'],'','Yes','No','No');
%         if ~strcmp(answer,'Yes'), disp('batch interrupted'), return, end
%     end
%     for k=kstart:S.nexp
%         fprintf('register trial %i/%i\n',k,S.nexp)
%         Y = fast_loaddata([S.datadir S.files(k,:)]);
%         S.ISO(:,:,k) = fast_recalage(Y,S.CS,iso);
%         iso = S.ISO(end,:,k);
%         save([S.savedir S.base '_batch'],'S')
%     end
%     
%     S.step = 'registered';
%     save([S.savedir S.base '_batch'],'S')
% end
% 
% 
% % 7 - speed estimation (preprocessing "on the fly")
% 
% if strcmp(S.step,'registered')
%     disp('speed estimation')
%     ns = length(S.edges);
% 
%     for kexp = 304:S.nexp
%         
%         timexp = tic;
%         fprintf('\n%s\n\nTRIAL %i/%i\n',datestr(now),kexp,S.nexp)
%         fname = sprintf('%s_%.2i',S.base,kexp);
%         
%         tic
%         if 1
%             disp('load')
%             Y = fast_loaddata([S.datadir S.files(kexp,:)]);
%             
%             Y = fast_recalage(Y,S.ISO(:,:,kexp));
%             %         fast_save2bytes(Y,[savedir fname '_register.ff']);
%         else
%             disp('load registered')
%             Y = fast_load2bytes([savedir fname '_register.ff']);
%         end
%         
%         Y = fast_normalizeseq(Y,[1 1 1 1 1]);
%         %     fast_save2bytes(Y,[savedir fname '_normalizeseq.ff']);
%         toc
%         
%         y = reshape(Y,S.nx*S.ny,S.nt);
%         
%         for k = 1:ns
%             tic
%             fprintf('vessel %i/%i\n',k,ns)
%             
%             I = S.edges(k).a * y;
%             I = I(:,6:end-5);
%             [J V] = fast_track(filty(I,20,'h'));
%             
%             e(kexp,k).I = I;
%             e(kexp,k).J = J;
%             e(kexp,k).V = V;
%             toc
%         end
%         
%         if mod(kexp,25)==0 || kexp==S.nexp
%             disp('save')
%             tic
%             save([S.savedir 'test_estimatedspeeds_part2'],'e','-V7.3')
%             toc
%         end
%         
%         fprintf('Trial elapsed time is %f seconds.',toc(timexp))
%     end
%     
% end



    