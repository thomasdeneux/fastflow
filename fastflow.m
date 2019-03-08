classdef fastflow < hgsetget
    properties
        nickname
    end
    properties (Transient)
        savedir = '';
    end
    properties
        datadir = '';
        resampledir = '';
        segmentpar
        parameters
        files
        filesoption = struct;
        cond
        nx
        ny
        nt
        CS
        nexp
        mask
        ISO
        edges = fast_edge.empty;
        flaglist = struct('k',1, ...
            'alllists',struct('name','flags list', ...
            'list',struct('flag','','color',[0 0 1],'defaultpars',[])));
        trials
        isdef = false;
    end
    properties (Access='private')
        version = 7; % removed the 'fake' field: now movie / line scan / fake edges can coexist!, added the 'isdef' field
        stim
        rest
    end
    properties (Transient, Access='private')
        fake % for backward compatibility only...
    end
    properties (Transient) %, Access='private')
        hf
        closefig
        ycontains % what does the global variable 'Y' contain: [#trial resampled?] or fakepar
        faketruespeed
    end
    properties (Transient)
        id_vol
    end
    
    % Constructor
    methods
        function S=fastflow
            % few initializations
            set(0,'defaultfigurecolormap',gray(64))
            fn_imvalue image
           
            % parameters
            S.parameters = fastflow.defaultparameters;
        end
    end
    
    % Buttons
    methods
        function closeproject(S)
            S.isdef = false;
            close(S.closefig(ishandle(S.closefig)))
        end
        function setoptions(S)
            [dum spec] = fastflow.defaultparameters;
            pars = fn_structedit(S.parameters,spec); 
            if isempty(pars), return, end
            S.parameters = pars;
            S.id_vol = [];
            saveproject(S)
        end
        function setedges(S)
            if ~status(S,'base')
                msgbox('base information missing')
                return
            end
            
            % put data into global array Y
            loadtrial(S,1);
            
            % edges selection:
            % - when 'SAVE' button will be pressed, the edges will be saved
            % - if a new project is opened, the figure will be closed
            E = fast_edgesmouse(S);
            S.closefig = union(S.closefig,E.hf);
        end
        function resultsvol(S)
            if ~status(S,'movie'), msgbox('you are not in movie mode'), return, end
            fvolname = [S.savedir 'volume' S.id_vol '.mat'];
            if ~exist(fvolname,'file'), msgbox('average volume time courses in image not computed yet'), return, end
            disp('load volume results'), drawnow
            v=load(fvolname);            
            if ~any(v.voldne), msgbox('average volume time courses in image not computed yet'), return, end
            %oi_blockview({volrest,volstim},1:5,[],[],[],[253 254],{'rest' 'stim'})
            figure(356), set(356,'numbertitle','off','name','RESULT VOLUME'), clf
            if S.parameters.resultavgtrials
                volall = (v.volstim+v.volrest)/2;
                volres = v.volstim./v.volrest;
                volres = fn_add(volres, ...
                    -mean(volres(:,:,1:S.parameters.volumefrzsub),3));
                avg = mean(volall,3);
                subplot(223), fn_4Dview(cat(4,fn_mult(v.volstim,1./avg),fn_mult(v.volrest,1./avg)),'2dplot'), title('stim+rest')
                subplot(224), fn_4Dview(volres,'2dplot'), title('stim/rest')
            else
                volall = v.volcat(:,:,:); 
                avg = mean(volall,3);
            end
            volres = fn_mult(volall,1./avg);
            clipavg = [min(avg(:)) max(avg(:))];
            clipres = [min(volres(:)) max(volres(:))];
            haim=subplot(121); fn_4Dview(avg,'2d','clip',clipavg)
            subplot(122), fn_4Dview(fn_mult(volall,1./avg),'2dplot'), title('all')
            hsens = fn_sensor('units','normalized','position',[.012 0.87 0.12 0.05], ...
                'value',get(haim,'clim'), ...
                'callback',@chgclip);
            hu = uicontrol('units','normalized','pos',[.012 0.93 0.12 0.05], ...
                'style','togglebutton','string','REF/RES', ...
                'callback',@togglerefres);
            function chgclip(u,evt) %#ok<INUSD>
                clip = get(hsens,'value');
                set(haim,'clim',clip)
                if get(hu,'value')
                    clipres = clip;
                else
                    clipavg = clip;
                end
            end
            function togglerefres(u,evt) %#ok<INUSD>
                if get(hu,'value')
                    data = volres;
                    clip = clipres;
                else
                    data = avg;
                    clip = clipavg;
                end
                hsens.value = clip;
                fn_4Dview('in',haim,'2d',data,'clip',clip)
            end
            S.closefig = union(S.closefig,356);
        end
        function results(S)
            % correction for possible mismatch btw. 'savedir' in S and in
            % S.edges
            if isempty(S.edges)
                msgbox('no eddges: define edges first and run analysis')
                return
            end
            R = fast_displayresult(S);
            S.closefig = union(S.closefig,R.hf);
        end
    end
    
    % Load/Save
    methods (Static)
        function [dpar spec] = defaultparameters
            dpar = struct('registration','isometry', ...
                'saveresampled',true, ...
                'volumexbin',2,'volumetbin',10,'volumefrzsub',10, ...
                'resultavgtrials',true);
            spec = struct('registration',{{'none' 'isometry(frame)'  'isometry(trial)' 'nonlinear(trial)'}}, ...
                'volumexbin','stepper 1 1 Inf 1', ...
                'volumetbin','stepper 1 1 Inf 1', ...
                'volumefrzsub','stepper 1 1 Inf 1');
        end
        function S = loadobj(S)
            % old versions
            if S.version<4
                msgbox('loading fastflow object from old version WILL make problem')
            end
            if S.version<5 
                if isfield(S.parameters,'fastmain')
                    S.parameters = S.parameters.fastmain;
                end
                % no more stim and rest, but condition number
                S.cond = zeros(1,S.nexp);
                S.cond(S.stim) = 1;
                S.stim = []; S.rest = [];
            end
            if S.version<5.1
                % indicate registration method rather than only true or
                % false
                doreg = double(S.parameters.doregistration);
                S.parameters = rmfield(S.parameters,'doregistration');
                S.parameters.registration = fn_switch(doreg, ...
                    0, 'none', 1, 'isometry(frame)', -1, 'isometry(trial)');
                if doreg==-1, S.ISO = S.ISO(1,:,:); end % registration values equal for all frames: keep only for first frame
            end
            if S.version<5.2
                S.fake = ~isempty(S.fake);
            end
            if S.version<5.3
                chgfield(S.edges,'nexp',S.nexp)
            end
            if S.version<5.4
                % now we can have several flag lists
                S.flaglist = struct('k',1,'alllists', ...
                    struct('name','flag list','list',S.flaglist));
            end
            if S.version<7
                S.isdef = true; % we asssumed that S is correctly defined!
            end
            S.version = 7;
            % use the current format for parameters
            S.parameters = fn_structmerge(fastflow.defaultparameters, ...
                S.parameters,'skip','recursive');
        end
    end
    methods
        function saveproject(S)
            x = S;
            % save the full version
            save([x.savedir x.nickname '_batch'],'x')
            % and a version without precomputed averages
            hidetmpdata(x.edges,'hide')
            save([x.savedir x.nickname '_light_batch'],'x')
            hidetmpdata(x.edges,'show')
        end
    end
    
    % Status
    methods
        function b = status(S,flag)
            switch flag
                case 'folders'
                    b = exist(S.savedir,'dir') && ~isempty(S.nickname);
                case 'datadir'
                    b = exist(S.datadir,'dir');
                case 'resampledir'
                    % if resampled data is not saved, status('resampledir')
                    % returns true even if the directory does not exist
                    b = ~S.parameters.saveresampled || exist(S.resampledir,'dir');
                case 'files'
                    b = status(S,'folders') && (~isempty(S.cond));
                case 'base'
                    b = status(S,'folders') && S.isdef;
                case 'movie'
                    b = ~isempty(S.mask) && ~isequal(S.mask,0); % in some old versions S.mask is set to 0 for fake
                case 'edges'
                    b = ~isempty(S.edges);
                otherwise
                    error('unknown status flag ''%s''',flag)
            end
        end
    end
    
    % Main routines
    methods
        function updatefolders(S)
            needdatadir   = status(S,'movie') || ~all(strcmp([S.edges.interpflag],'fake'));
            needresampdir = needdatadir && S.parameters.saveresampled && any(strcmp({S.edges.interpflag},'movie'));
            if needdatadir && ~exist(S.datadir,'dir')
                S.datadir = fn_getfile('DIR','Select data directory'); 
                if isequal(S.datadir,0), S.datadir=''; else S.datadir(end+1)='/'; end
            end
            if needresampdir && ~exist(S.resampledir,'dir')
                S.resampledir = fn_getfile('DIR','Select resampled data directory'); 
                if isequal(S.resampledir,0), S.resampledir=''; else S.resampledir(end+1)='/'; end
            end
        end
        function inputwidefield(S)
            global Y
            % 1 - data and saving directory
            if ~status(S,'folders')
                disp('Select data directory')
                S.datadir = fn_getfile('DIR','Select data directory');
                if isequal(S.datadir,0), disp('batch interrupted'), return, end
                S.datadir(end+1) = '/';
                disp('Select save directory')
                S.savedir = fn_getfile('DIR','Select save directory');
                if isequal(S.savedir,0), disp('batch interrupted'), return, end
                S.savedir(end+1) = '/';
                disp('Select resampled data directory')
                S.resampledir = fn_getfile('DIR','Select resampled data directory');
                if isequal(S.resampledir,0), disp('batch interrupted'), return, end
                S.resampledir(end+1) = '/';
                
                if isempty(S.nickname)
                    [dum default] = fileparts(S.datadir); %#ok<*ASGLU>
                else
                    default = S.nickname;
                end
                answer = inputdlg('Enter base nickname','',1,{default});
                if isempty(answer), disp('batch interrupted'), return, end
                S.nickname = answer{1};
                %saveproject(S)
            end
            if ~status(S,'datadir'), disp('could not find data directory'), end
            if ~status(S,'resampledir'), disp('could not find resampled data directory'), end
            
            % 2 - data files, conditions, sort according to experiments and conditions
            if ~status(S,'files')
                % data files
                disp('file names')
                b = dir(S.datadir);
                a = [];
                f = false(1,length(b));
                
                % try format BLABLAC01_E01B001.BLK
                if isempty(a)
                    for k=1:length(b)
                        tokens = regexp(b(k).name,'.*C(\d{2}).*_E(\d{2})B(\d{2,3}).BLK$','tokens');
                        if isempty(tokens), continue, end
                        f(k) = true;
                        b(k).cond  = str2double(tokens{1}{1});
                        b(k).exp   = str2double(tokens{1}{2});
                        b(k).block = str2double(tokens{1}{3});
                    end
                    a = b(f);
                end
                
                % try format BLABLA_E01B001.BLK
                if isempty(a)
                    for k=1:length(b)
                        tokens = regexp(b(k).name,'.*_E(\d{2})B(\d{2,3}).BLK$','tokens');
                        if isempty(tokens), continue, end
                        f(k) = true;
                        b(k).cond  = 0;
                        b(k).exp   = str2double(tokens{1}{1});
                        b(k).block = str2double(tokens{1}{2});
                    end
                    a = b(f);
                end
                
                % try format BLABLA01.a
                if isempty(a)
                    nstim = 0;
                    for k=1:length(b)
                        tokens = regexp(b(k).name,'.*(\d{2})\.([a-z])$','tokens');
                        if isempty(tokens), continue, end
                        if ~nstim % first time
                            nstim = get_head([S.datadir b(k).name],'nstimuli');
                        end
                        f(k) = true;
                        b(k).exp   = str2double(tokens{1}{1});
                        b(k).block = str2double(tokens{1}{2})-(double('a')-1);
                    end
                    f = find(f);
                    a = b(f)';
                    a = repmat(a,nstim,1);
                    for i=1:nstim
                        [a(i,:).cond] = deal(i); 
                        condstr = num2str(i);
                        for k=1:length(f), a(i,k).name = [b(f(k)).name '-cond' condstr]; end
                    end
                    a = a(:)';
                end
                
                % try mat files
                if isempty(a)
                    a = dir([S.datadir 'trial*_cond*.mat']);
                    for k=1:length(a)
                        tokens = regexp(a(k).name,'trial(\d*)_cond(\d*)','tokens');
                        tokens = tokens{1};
                        a(k).exp   = 1;
                        a(k).block = k;
                        a(k).cond  = str2double(tokens{2});
                    end
                end
                
                % try avi files
                if isempty(a)
                    a = dir([S.datadir '*.avi']);
                    if ~isempty(a)
                        [a.exp] = deal(1);
                        for k=1:length(a), a(k).block = k; end
                        [a.cond] = deal(0);
                    end
                end
                
                if isempty(a)
                    error('cannot find data files: try defining another filter')
                end
                
                % select which experiments to take
                exps = unique([a.exp]);
                answer = inputdlg('Select experiments','',1,{num2str(exps)});
                if isempty(answer), disp('batch interrupted'), return, end
                exps = str2num(answer{1}); %#ok<ST2NM>
                f = ismember([a.exp],exps);
                a = a(f);
                % sort files
                exps   = [a.exp]';
                blocks = [a.block]';
                [dum ord] = sortrows([exps blocks]);
                a = a(ord);
                hfg = figure('defaultuicontrolunits','normalized'); clf
                uicontrol('pos',[.02 .12 .96 .78],'style','listbox','string',{a.name});
                uicontrol('pos',[.02 .02 .6 .07],'style','text','string','Approve files order?');
                h1=uicontrol('pos',[.64 .02 .15 .07],'string','OK','callback',@(u,evnt)delete(u));
                h2=uicontrol('pos',[.83 .02 .15 .07],'string','Cancel','callback',@(u,evnt)delete([h1 u]));
                waitfor(h1)
                ok = ishandle(h2); % if OK pressed, h1 is deleted but not h2
                if ishandle(hfg), close(hfg), end
                if ~ok, disp('batch interrupted'), return, end
                S.files = strvcat(a.name); %#ok<VCAT>
                S.nexp = size(S.files,1);
                
                % check for color video
                fname = deblank(S.files(1,:)); % first file
                if strcmpi(fn_fileparts(fname,'ext'),'.avi')
                    v = VideoReader(fullfile(S.datadir,fname));
                    if ~strfind(v.VideoFormat,'RGB')
                        error('Currently only video format with 3 color channels is handled. Please edit code for new format.')
                    end
                    answer = questdlg('Video appear to have 3 color channels. How do you want to convert it to grayscale?', ...
                        '', 'Take the first channel', 'Average over 3 channels (more memory needed)', 'Take the first channel');
                    if isempty(answer)
                        disp 'interrupted'
                        return
                    end
                    S.filesoption.average_rgb = strcmp(answer, 'Average over 3 channels (more memory needed)');    
                end
                
                % conditions
                conds = [a.cond]';
                tmp = unique(conds)';
                answer = inputdlg({'rest conditions:','stim1 conditions:','stim2 conditions:'}, ...
                    'Define conditions', ...
                    1, ...
                    {num2str(tmp([1 end])),num2str(tmp(2:end-1)),''});
                condstim1 = str2num(answer{2}); %#ok<ST2NM>
                stimtrials1 = ismember(conds,condstim1);
                condstim2 = str2num(answer{3}); %#ok<ST2NM>
                stimtrials2 = ismember(conds,condstim2);
                S.cond = zeros(1,S.nexp); 
                S.cond(stimtrials1) = 1;
                S.cond(stimtrials2) = 2;
                %saveproject(S)
            end
            
            % 3 - basic data
            if ~status(S,'base')               
                loadrawdata(S,1)
                %                 fname = deblank([S.datadir S.files(1,:)]);
                %                 [dum1 dum2 ext] = fileparts(fname); %#ok<ASGLU>
                %                 switch ext(2:end)
                %                     case 'BLK'
                %                         fname = regexprep(fname,'-cond\d$','');
                %                         hdr = fast_loaddata(fname,'header');
                %                         S.nx = hdr.xs;
                %                         S.ny = hdr.ys;
                %                         S.nt = hdr.nfrms;
                %                     case 'mat'
                %                         v = load(fname); F = fieldnames(v);
                %                         if ~isscalar(F), error('wrong data file'), end
                %                         data = v.(F{1}); if ndims(data)~=3, error('wrong data'), end
                %                         [S.nx S.ny S.nt] = size(data);
                %                     otherwise
                %                         error('cannot open data file with extension ''%s''',ext(2:end))
                %                 end
                %                 fast_loaddata_global(fname)
                [S.nx S.ny S.nt] = size(Y);
                CS_ = mean(Y,3);
                answer = '';
                while ~strcmp(answer,'Yes')
                    mask_ = fast_structure(CS_);
                    figure(1), imagesc(CS_')
                    figure(2), imagesc(mask_')
                    answer = questdlg('Approve reference frame?','','Yes','No','Quit','No');
                    close(1:2)
                    if strcmp(answer,'Quit')
                        return
                    elseif strcmp(answer,'No')
                        figure(1), fn_4Dview(Y,'2d')
                        figure(2), fn_4Dview(Y,'2dplot')
                        frames = fn_input('Select frames to average for reference',1:10);
                        close(1:2)
                        CS_ = mean(Y(:,:,frames),3);
                    end
                end
                S.CS = CS_;
                S.mask = mask_;
                % selected trials
                S.trials = true(1,S.nexp);
                % registration parameters
                S.ISO = zeros(S.nt,1,S.nexp);
                S.isdef = true;
                saveproject(S)
            end
        end
        function inputlinescan(S)
            % 1 - saving directory
            if ~status(S,'folders')
                disp('Select data directory')
                S.datadir = fn_getfile('DIR','Select data directory');
                if isequal(S.datadir,0), disp('batch interrupted'), return, end
                S.datadir(end+1) = '/';
                disp('Select save directory')
                S.savedir = fn_getfile('DIR','Select save directory');
                if isequal(S.savedir,0), S.savedir = ''; disp('batch interrupted'), return, end
                S.savedir(end+1) = '/';
                
                [dum base] = fileparts(S.savedir(1:end-1));
                answer = inputdlg('Enter base nickname','',1,{base});
                if isempty(answer), disp('batch interrupted'), return, end
                S.nickname = answer{1};
                saveproject(S)
            end
            if ~status(S,'datadir'), disp('could not find data directory'), end

            % 2 - data files, conditions, sort according to experiments and conditions
            if ~status(S,'files')
                % data files
                disp('file names')
                b = dir(S.datadir);
                a = [];
                f = false(1,length(b));
                
                % try format
                % linescan-015_Cycle00001_CurrentSettings_Ch2_000198.tif
                if isempty(a)
                    for k=1:length(b)
                        tokens = regexp(b(k).name,'^linescan-(\d)*_Cycle\d*_CurrentSettings_Ch\d_(\d)*\.tif$','tokens');
                        if isempty(tokens), continue, end
                        f(k) = true;
                        b(k).cond  = 0;
                        b(k).edge  = str2double(tokens{1}{1});
                        b(k).block = str2double(tokens{1}{2});
                    end
                    a = b(f);
                end
                
                if isempty(a)
                    error('cannot find data files: try defining another filter')
                end
                
                % select which edges to take
                edgesq = unique([a.edge]);
                answer = inputdlg('Select scans','',1,{num2str(edgesq)});
                if isempty(answer), disp('batch interrupted'), return, end
                edgesq = str2num(answer{1}); %#ok<ST2NM>
                f = ismember([a.edge],edgesq);
                a = a(f);
                % sort files
                edges  = [a.edge]';
                blocks = [a.block]';
                [dum ord] = sortrows([edges blocks]);
                a = a(ord);
                hfg = figure('defaultuicontrolunits','normalized'); clf
                uicontrol('pos',[.02 .12 .96 .78],'style','listbox','string',{a.name});
                uicontrol('pos',[.02 .02 .6 .07],'style','text','string','Approve files order?');
                h1=uicontrol('pos',[.64 .02 .15 .07],'string','OK','callback',@(u,evnt)delete(u));
                h2=uicontrol('pos',[.83 .02 .15 .07],'string','Cancel','callback',@(u,evnt)delete([h1 u]));
                waitfor(h1)
                ok = ishandle(h2); % if OK pressed, h1 is deleted but not h2
                if ishandle(hfg), close(hfg), end
                if ~ok, disp('batch interrupted'), return, end
                
                % create edges
                ne = length(edgesq);
                for k=1:ne
                    ak = a;
                    S.edges(k) = fast_edge(S,'scan',strvcat(ak.name)); %#ok<*VCAT>
                end
                
                % save
                S.isdef = true;
                saveproject(S)
            end
        end
        function inputfake(S)
            % 1 - saving directory
            if ~status(S,'folders')
                S.datadir = [];
                disp('Select save directory')
                S.savedir = fn_getfile('DIR','Select save directory');
                if isequal(S.savedir,0), S.savedir = ''; disp('batch interrupted'), return, end
                S.savedir(end+1) = '/';
                
                [dum base] = fileparts(S.savedir(1:end-1));
                answer = inputdlg('Enter base nickname','',1,{base});
                if isempty(answer), disp('batch interrupted'), return, end
                S.nickname = answer{1};
                saveproject(S)
            end
            
            % 2 - fake parameters
            if ~status(S,'base')
                %  basic data
                S.files = [];
                S.nexp = 0;
                S.cond = [];
                S.rest = [];
                % selected trials
                S.trials = true(1,S.nexp);
                % registration parameters
                S.parameters.registration = 'none';
                S.isdef = true;
                saveproject(S)
            end
        end
        function analysis(S,flag)
            if nargin<2, flag='all'; end
            if ~fn_ismemberstr(flag,{'all','interp+vol'})
                error('unknown analysis flag ''%s''',flag)
            end
            
            % global array
            global Y 
            
            % local parameters
            pars = S.parameters;
                
            % few checks
            if ~status(S,'base')
                msgbox('base information missing')
                return                
            end
            if all(size(S.edges)>1)
                error('edges must be a vector')
            end
            
            % create/load file of saved resampled volume
            if status(S,'movie')
                disp('check volume file'), drawnow
                fvolname = [S.savedir 'volume' S.id_vol '.mat'];
                okfile = exist(fvolname,'file');
                if ~okfile
                    voldne = false(1,S.nexp); 
                    volbin = [pars.volumexbin pars.volumetbin];
                    siz = floor([S.nx S.ny S.nt]./volbin([1 1 2]));
                    disp('new volume data file')
                    if pars.resultavgtrials
                        volstim = zeros(siz,'single');
                        volrest = zeros(siz,'single');
                        save(fvolname,'volstim','volrest','volbin','voldne')
                    else
                        siz(4)=S.nexp;
                        volcat = zeros(siz,'single');
                        save(fvolname,'volcat','volbin','voldne')
                    end
                else
                    load(fvolname,'voldne')
                    if ~all(voldne) %#ok<NODEF>
                        load(fvolname)
                    end
                end
            end
            
            % GO! - loop on trials
            activeedge = [S.edges.active];
            ntrial = max([S.nexp S.edges(activeedge).nexp]); % in fake mode, nexp can differ between edges
            for ktrial=1:ntrial
                fprintf('\nTRIAL %i/%i\n',ktrial,ntrial)
                tim = tic;
                
                % 1 - coregistration
                if status(S,'movie') && ~strcmp(pars.registration,'none') ...
                        && status(S,'datadir') && status(S,'resampledir')
                    fname = [S.resampledir strrep(deblank(S.files(ktrial,:)),'.BLK','') '.ff'];
                    fname = deblank(fname);
                    if ~any(any(S.ISO(:,:,ktrial))) ...
                            || (pars.saveresampled && ~exist(fname,'file'))
                        % registration not done yet, or registration done
                        % but still need to resample the data
                        % the routine loadtrial(S,ktrial) in fact does all
                        % what is needed
                        loadtrial(S,ktrial)
                    end
                end                
                
                % 2 - resample the whole movie (volume)
                if status(S,'movie') 
                    if ~voldne(ktrial)
                        loadtrial(S,ktrial);
                        disp('bin data')
                        vol = fn_bin(Y,volbin([1 1 2]));
                        disp('save binned')
                        voldne(ktrial) = true;
                        if S.parameters.resultavgtrials
                            if S.cond(ktrial)>0
                                volstim = volstim + vol/sum(S.cond>0);
                            else
                                volrest = volrest + vol/sum(S.cond==0);
                            end
                            save(fvolname,'volstim','volrest','volbin','voldne')
                        else
                            volcat(:,:,:,ktrial) = vol;
                            save(fvolname,'volcat','volbin','voldne')
                        end
                    end
                end
                
                % 3 - interpolate vessel data and vessel section data
                okedge = [S.edges.active] & ([S.edges.nexp]>=ktrial); % in fake mode, nexp can differ between edges
                idxokedge = find(okedge);
                edgesk = S.edges(okedge);
                bedge = hasdata(edgesk(:),ktrial)';
                sectionsk = [S.edges([S.edges.nexp]>=ktrial).section];
                sectionsk = sectionsk([sectionsk.active]);
                bsection = hassectiondata(sectionsk,ktrial)';
                if any(~[bedge bsection])
                    if status(S,'movie'), loadtrial(S,ktrial); end
                    disp(['vessels data (' num2str(sum(~[bedge bsection])) ' vessels and vessel sections)'])
                    drawnow
                    for i=find(~bedge)
                        e = edgesk(i);
                        switch e.interpflag
                            case 'movie'
                                setdata(e,ktrial,fast_interptube(Y,e));
                            case 'scan'
                                fname = [S.datadir deblank(e.filepar.files(ktrial,:))];
                                setdata(e,ktrial,fn_readimg(fname));
                            case 'fake'
                                setdata(e,ktrial,fast_fakedata(e.fakepar));
                        end
                    end
                    for i=find(~bsection)
                        ui = sectionsk(i);
                        setsectiondata(ui,ktrial,fast_vesselsection(Y,ui.edge,ui.sectionpar));
                    end
                end
                
                % 4 - estimate flow (only on vessels selected for analysis)
                if strcmp(flag,'all')
                    ne = numel(edgesk);
                    for i=1:ne
                        e  = edgesk(i);
                        if ~e.active, error programming, end
                        nf = numel(e.flow);
                        firsttime = true;
                        for j=1:nf
                            flow = e.flow(j);
                            if ~flow.active || hasresult(flow,ktrial), continue, end
                            % first time action
                            if firsttime
                                % display
                                fprintf('\nestimate flow, vessel %i (%i/%i)\n',idxokedge(i),i,ne)
                                %if ~isempty(e.fakepar), disp(e.fakepar), end
                                firsttime = false;
                            end
                            % data: filtering
                            I = getdataf(flow,ktrial);
                            % estimate
                            switch flow.method
                                case 'track'
                                    par = flow.estpar;
                                    % check special initialization
                                    J0V0={};
                                    if ischar(par.v0)
                                        switch par.v0
                                            case 'last'
                                                if ktrial==1
                                                    par.v0 = [];
                                                else
                                                    V0 = getresult(flow,ktrial-1);
                                                    par.v0 = mean(V0(:));
                                                end
                                            case 'radon'
                                                % first get a radon estimate!
                                                methods = {e.flow.method}; nps = [e.flow.np]; nts = [e.flow.nt];
                                                kradon = find(strcmp(methods,'radon') & nps==flow.np & nts==flow.nt,1);
                                                if isempty(kradon)
                                                    disp('no radon result to initialize track')
                                                    continue
                                                end
                                                fradon = e.flow(kradon);
                                                if ~hasresult(fradon,ktrial)
                                                    disp(fradon.algo), tic
                                                    V0 = fast_radon(I,fradon.estpar);
                                                    setresult(fradon,ktrial,V0)
                                                    fradon.user.time(ktrial) = toc;
                                                    fprintf('\b (%is)\n',floor(toc))
                                                end
                                                V0 = getresult(fradon,ktrial);
                                                J0V0 = {I,V0};
                                            case {'fft', 'est'}
                                                % keep it
                                            otherwise
                                                error merde
                                        end
                                    end
                                    disp(flow.algo), tic
                                    [J V] = fast_track(I,par,J0V0{:});
                                    setresult(flow,ktrial,V,J)
                                    flow.user.time(ktrial) = toc;
                                    fprintf('\b (%is)\n',floor(toc))
                                case 'gabor'
                                    disp(flow.algo), tic
                                    V = fast_gabor(I,flow.estpar);
                                    setresult(flow,ktrial,V)
                                    flow.user.time(ktrial) = toc;
                                    fprintf('\b (%is)\n',floor(toc))
                                case 'radon'
                                    disp(flow.algo), tic
                                    V = fast_radon(I,flow.estpar);
                                    setresult(flow,ktrial,V)
                                    flow.user.time(ktrial) = toc;
                                    fprintf('\b (%is)\n',floor(toc))
                            end
                        end
                    end
                end
                
                % 5 - estimate vessel width (only on vessels selected for analysis)
                if strcmp(flag,'all')
                    b = hasresult(sectionsk,ktrial)';
                    fb = find(~b);
                    nb = length(fb);
                    fprintf('\n')
                    for i=1:nb
                        j = fb(i);
                        ui = sectionsk(j);
                        % display
                        fprintf('estimate section width %i (%i/%i)\n',j,i,nb)
                        % estimate
                        if ktrial==1, p0={}; else p0={mean(getresult(ui,struct('data','profile','trial',ktrial-1)),2)}; end
                        profile = fast_vesselwidth(getsectiondata(ui,ktrial),p0{:});
                        setresult(ui,ktrial,profile);
                    end
                end
                
                % display
                tim = toc(tim);
                fprintf('\nElapsed time for TRIAL %i: ',ktrial)
                if tim<60
                    fprintf('%.1fs\n\n',tim)
                else
                    fprintf('%.1fmin\n\n',tim/60)
                end
                remains = (ntrial-ktrial)*tim;  % seconds
                pred = now + remains/(24*3600); % days          
                if floor(pred)-floor(now)>1
                    fprintf('\bEnd expected %s\n\n',datestr(pred,'dd mmm, HH:MM'))
                elseif floor(pred)-floor(now)==1
                    fprintf('\bEnd expected tomorrow, %s\n\n',datestr(pred,'HH:MM'))
                elseif remains>660
                    fprintf('\bEnd expected at %s\n\n',datestr(pred,'HH:MM'))
                elseif remains>60
                    fprintf('\bEnd expected in %i min\n\n',floor(remains/60))
                end
            end
            
            % Average results over trials
            if strcmp(flag,'all')
                precomp(S)
            end
        end
        function precomp(S,ktrial)
            % input
            if nargin<2, ktrial=[]; end
            if length(ktrial)>1, error('ktrial must be scalar'), end
            % prepare resdisp structure
            resdisp = struct('data',[],'algo',[], ...
                'trial',[],'res',[],'oddeven',[], ...
                'points',[],'fzsub',0,'smooth',0);
            ne = numel(S.edges);
            for ke = 1:ne
                e = S.edges(ke);
                %if ~e.active, continue, end
                fprintf('precompute average for vessel %i/%i\n',ke,ne)
                nf = length(e.flow);
                for i=0:nf
                    % volume or result?
                    if i==0
                        obj = e;
                        resdisp.data = 'volume';
                        resdisp.algo = 'volume';
                        datas = {'volume'};
                        hasdata(obj,[]);
                    else
                        obj = e.flow(i);
                        resdisp.data = 'result';
                        resdisp.algo = '';
                        datas = {'result' 'lines'};
                        hasresult(obj,[]);
                    end
                    if ~obj.active, continue, end
                    % average results
                    if S.parameters.resultavgtrials
                        resdisp.trial   = 'average';
                        resdisp.res     = 'all';
                        resdisp.oddeven = 'all';
                        getresult(obj,resdisp);
                    else
                        resdisp.trial = 'timeavg';
                        getresult(obj,resdisp);
                    end
                    % specific trial results
                    if ~isempty(ktrial)
                        resdisp.trial = ktrial;
                        for data = datas
                            resdisp.data = data{1};
                            getresult(obj,resdisp);
                        end
                    end
                end
            end
        end
    end
    
    % Results
    methods
        function id = get.id_vol(S)
            par = S.parameters;
            s = struct('FILES',{S.files}, ...
                'XBIN',par.volumexbin,'TBIN',par.volumetbin, ...
                'AVGTRIALS',par.resultavgtrials);
            S.id_vol = fn_hash(s,4);
            fil = ['volume' S.id_vol '.xml'];
            if exist(S.savedir,'dir') && ~exist([S.savedir fil],'file')
                disp(['write estimation parameters in file ' fil])
                fn_savexml([S.savedir fil],s);
            end
            id = S.id_vol;
        end
        function set.trials(S,trials)
            S.trials = trials;
        end
        function x = getresult(S,resdisp)
            x = getresult(S.edges([S.edges.active]),resdisp);
        end
    end
    
    % Flag lists
    methods
        function [list name] = getflaglist(S)
            list = S.flaglist.alllists(S.flaglist.k);
            if nargout==2, name = list.name; end
            list = list.list;
        end
        function setflaglist(S,varargin)
            % function setflaglist(S,list)
            % function setflaglist(S,k,s)
            klist = S.flaglist.k;
            switch nargin
                case 2
                    list = varargin{1};
                    S.flaglist.alllists(klist).list = list;
                case 3
                    [k s] = deal(varargin{:});
                    S.flaglist.alllists(klist).list(k) = s;
                otherwise
                    error argument
            end
        end
    end
    
    % Files control
    methods
        function deletefiles = cleanfiles(S)
            swd = pwd;
            cd(S.savedir)
            % find which files represent valid data
            d = [dir('edge*_*'); dir('fake*_*'); dir('scan*_*')];
            d = d(~[d.isdir]);
            allfiles = {d.name};
            goodfiles = {[S.nickname '_batch.mat'], [S.nickname '_light_batch.mat']};
            if ~isempty(S.edges)
                allflows = [S.edges.flow];
                allsections = [S.edges.section];
                goodfiles = [goodfiles S.edges.files allflows.files allsections.files];
            end
            okfile = ismember(allfiles,goodfiles);
            deletefiles = allfiles(~okfile);
            % return the list or delete files
            if nargout==0
                if isempty(deletefiles)
                    disp('NO FILE TO DELETE')
                else
                    disp('FILES TO DELETE')
                    nfile = length(deletefiles);
                    for i=1:nfile, disp(deletefiles{i}), end
                    answer = input('do really delete files (y/n)? ','s');
                    if any(strcmp(answer,{'y' 'yes'}))
                        delete(deletefiles{:})
                    end
                end
                clear deletefiles
            end
            cd(swd)
        end
    end
    
    % User tools
    methods
        %         function createfake(S,s)
        %             % function createfake(S,fakepar)
        %             global Y
        %             [Y x] = fast_fakedata(s);
        %             m = ~isnan(x);
        %             sigma = 5;
        %             x(~m) = 0;
        %             x = filt2(x,sigma)./filt2(m,sigma);
        %             S.faketruespeed = x; % TODO: handle true speed...
        %         end
        function success = loadrawdata(S,ktrial)
            global Y
            % function loadrawdata(S,ktrial)
            success = false;
            if isequal(S.ycontains,[ktrial 0])
                % global variable 'Y' already contains the desired data,
                % nothing to do
                success = true;
                return
            end
            S.ycontains = []; % in case function will be interrupted in the middle
            fname = [S.datadir deblank(S.files(ktrial,:))];
            % how to read file
            [dum1 dum2 ext] = fileparts(fname);  %#ok<ASGLU>
            ftype = ext(2:end);
            if strcmp(ftype,'BLK')
                token = regexp(fname,'-(cond\d)','tokens');
                if ~isempty(token)
                    % each file contains several trials (several
                    % conditions)
                    ftype = 'BLKcond';
                    fname = regexprep(fname,'-cond\d','');
                end
            end
            if ~exist(fname,'file'),
                if nargout==1,return, else error('cannot find raw data file'), end
            end
            switch lower(ftype)
                case 'blk'
                    fast_loaddata_global(fname);
                case 'blkcond'
                    Y = oi_loadBLK(fname,'array','double',token{1}{1});
                case 'mat'
                    v = load(fname);
                    F = fieldnames(v);
                    if ~isscalar(F), error('wrong data file'), end
                    Y = v.(F{1});
                    if ~strcmp(class(Y),'double'), Y = single(Y); end
                case 'avi'
                    x = fn_readmovie(fname, 'nopermute');
                    if S.filesoption.average_rgb
                        x = squeeze(mean(x,3)); % transform color movie into unidim movie
                    else
                        x = squeeze(x(:,:,1,:));
                    end
                    if isa(x,'double'), x = single(x); end
                    Y = permute(x, [2 1 3]);
                otherwise
                    error('cannot open file of type ''%s''',ftype)
            end
            % finish
            S.ycontains = [ktrial 0];
            if nargout==0
                clear success
            else
                success = true;
            end
        end
        function success = loadtrial(S,ktrial)
            % function loadtrial(S,ktrial)
            global Y
            success = false;
            if ~status(S,'movie')
                if ~success && nargout==0, error 'you are not in movie mode', end
                return
            end
            %ktrial = argh;
            regmethod = S.parameters.registration;
            doreg = ~strcmp(regmethod,'none');
            containflag = [ktrial doreg];
            if isequal(S.ycontains,containflag)
                % global variable 'Y' already contains the desired data,
                % nothing to do
                success = true;
                return
            end
            if ~doreg
                % no need for registration
                success = loadrawdata(S,ktrial); 
                if ~success && nargout==0, error('cannot find raw data file'), end
                return
            else
                if ~status(S,'resampledir'), disp('could not find resampled data directory'), end
                fname = [S.resampledir strrep(deblank(S.files(ktrial,:)),'.BLK','') '.ff'];
                if exist(fname,'file')
                    % resampled data exists - load it
                    disp('load resampled data'), drawnow
                    S.ycontains = []; % if computation will be interrupted, data in global Y might be corrupted
                    fast_load2bytes_global(fname);
                    S.ycontains = containflag;
                else
                    success = loadrawdata(S,ktrial);
                    if ~success
                        if nargout==0, error('could not find raw data file'), else return, end
                    end
                    if ~any(any(S.ISO(:,:,ktrial)))
                        % registration not done - do it
                        if ktrial==1
                            iso0 = [];
                        else
                            if findstr(regmethod,'isometry')
                                iso0 = S.ISO(end,:,ktrial-1);
                            elseif findstr(regmethod,'nonlinear')
                                iso0 = S.ISO(:,:,ktrial-1);
                            else
                                error('wrong registration method ''%s''',regmethod)
                            end
                            if ~any(iso0(:))
                                if nargout==0
                                    error('attempt to coregister trial %i before trial %i',ktrial,ktrial-1)
                                else
                                    return
                                end
                            end
                        end
                        S.ycontains = []; % in computation will be interrupted, data in global Y might be corrupted
                        x = fast_recalage_global(S.CS,regmethod,iso0);
                        if ktrial==1
                            % might change the size of S.ISO
                            [m n] = size(x);
                            S.ISO = zeros(m,n,S.nexp);
                            % avoid zero
                            x = x+1e-10;
                        end
                        S.ISO(:,:,ktrial) = x(:,:);
                        saveproject(S)
                    else
                        % registration done but not resampling - do it               
                        S.ycontains = []; % in computation will be interrupted, data in global Y might be corrupted
                        fast_recalage_global(regmethod,S.ISO(:,:,ktrial));
                    end
                    S.ycontains = containflag;
                    if S.parameters.saveresampled && status(S,'resampledir')
                        disp('save resampled data'), drawnow
                        fast_save2bytes(Y,fname)
                    end
                end
            end
            if nargout==0
                clear success
            else
                success = true;
            end
        end
        function Z = gettrial(S,ktrial,varargin)
            % function Z = gettrial(S,ktrial[,frames][,classname])
            global Y
            frames = 1:S.nt; classname = 'double';
            for i=1:length(varargin)
                a = varargin{i};
                if isnumeric(a), frames = a; else classname = a; end
            end
            if ~status(S,'movie'), error('you are not in movie mode'), end
            fname = [S.resampledir strrep(S.files(ktrial,:),'.BLK','.ff')];
            if ~S.parameters.doregistration || ~S.parameters.saveresampled 
                error('this part is not implemented yet')
            end
            % create resampled data if does not exist
            if ~exist(fname,'file')
                loadtrial(S,ktrial); 
                Z = zeros(S.nx,S.ny,length(frames),classname);
                for k=1:length(frames), Z(:,:,k) = Y(:,:,frames(k)); end
            else
                Z = fast_load2bytes(fname,frames,classname);
            end
        end
        function access(S) %#ok<MANU>
            keyboard
        end
    end
end





