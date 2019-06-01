classdef madflow < interface
   
    properties (SetAccess='private')
        S
    end
    % current diplay
    properties (SetAccess='private')
        ktrial
        kedge
        e = fast_edge.empty;
        kflow 
        f = fast_flow.empty;
        ksection
        u = fast_section.empty;
    end
    % structures with specific information
    properties (SetAccess='private')
        misc        % special information 
        im          % vessels image (wide field mode)
        vsl         % vessels (wide field mode)
        scan        % list of vessels (line scan and fake modes)
        polysel     % new vessel selection
        parest      % estimation parameters
        sections    % vessel section parameters
        trials      % trials
        res         % results (one element per graph)
        resd        % results (others)
    end
    
    % Constructor
    methods
        function M = madflow(S)
            % init interface
            defaultoptions = struct( ...
                'segment',      {fast_segment('parameters')}, ...
                'filterpar',    {{'y = y ./ filtx(y,20,''lm'');'
                                  'y = filty(y,10,''hm'');'
                                  'y = filty(y,1,''lm'');'}}, ...
                'parest',       struct, ...
                'sections',     struct, ...
                'resdisp',      repmat(fast_resdispselect('default'),1,7));
            M = M@interface(862,'FASTFLOW',defaultoptions);
            set(M.hf,'defaultuicontrolfontsize',10)
            
            % graphic objects
            init_grob(M)
            interface_end(M)
            drawnow
           
            % specific initializations
            misc_init(M)
            vessels_init(M)
            parest_init(M)
            sections_init(M)
            trials_init(M)
            res_init(M)
            
            % content
            if nargin<1
                S = fastflow;
            end
            M.S = S;
            fastflow_manage(M,'start')
        end
        function init_grob(M)
            % GENERAL
            bgcolor = [.8 .8 .7];
            M.grob.panel = axes('color',bgcolor,'xcolor',bgcolor,'ycolor',bgcolor, ...
                'xtick',[],'ytick',[],'hittest','off');
            info = {'vessels','parameters','trials','results','sections'};
            for i=1:length(info)
                M.grob.info(i) = uicontrol('style','text','string',info{i});
            end
            M.grob.actionbutton = uicontrol('backgroundcolor',bgcolor);
            
            % VESSELS
            % (image)
            M.grob.haim = axes('xtick',[],'ytick',[],'box','on');
            M.grob.im_toggle = togglegroup({'vessels','energy','movie'},@(s)im_toggle(M));
            M.im.objects = [M.grob.haim M.grob.im_toggle];
            % (list of vessels)
            M.grob.scanlabel = uicontrol('style','text','string','(list)');
            M.grob.scanlist = uicontrol('style','listbox','max',2, ...
                'callback',@(u,e)scan_manage(M,'select',get(u,'value')));            
            M.grob.scanorder  = uicontrol('string','reorder', ...
                'callback',@(u,evt)scan_manage(M,'reorder'));
            M.grob.scannewfake  = uicontrol('string','new fake', ...
                'callback',@(u,evt)scan_manage(M,'newfake'));
            M.scan.objects = [M.grob.scanlabel M.grob.scanlist ...
                M.grob.scanorder M.grob.scannewfake]; 
            % (vessels)
            M.misc.huvessel = fn_stepper('min',0,'max',0,'step',1, ...
                'format','%i','coerce',true,'value',0, ...
                'callback',@(u,evt)vessels_setcurrent(M,get(u,'value')));
            M.grob.huvessel = M.misc.huvessel.hpanel; % need the grob value to be a Matlab graphic object, not a fn_stepper object
            
            % PARAMETERS
            % (flags & colors)
            M.grob.haflags = axes('xtick',[],'ytick',[],'box','on');
            % (manage parameters list)
            M.grob.parlist = uicontrol('style','listbox', ...
                'callback',@(u,evt)parest_manage(M,'select',get(u,'value')));
            M.grob.parlisttop  = uicontrol('string','move to top', ...
                'callback',@(u,evt)parest_manage(M,'top'));
            M.grob.parlistrm   = uicontrol('string',' remove ', ... MATLAB BUG! 'remove' produces blank display! use ' remove ' instead
                'callback',@(u,evt)parest_manage(M,'remove'));
            M.grob.parlistsave = uicontrol('string','add', ...
                'callback',@(u,evt)parest_manage(M,'add'));
            M.grob.parlistreplace = uicontrol('string','replace', ...
                'callback',@(u,evt)parest_manage(M,'replace'));
            % (define parameters)
            M.grob.parmethod = uicontrol('style','popupmenu','string','x','value',1);
            M.grob.parid = uicontrol('style','text');
            M.grob.parfilter = uicontrol('style','edit','max',2, ...
                'backgroundcolor','white','horizontalalignment','left', ...
                'callback',@(u,evt)parest_manage(M,'changefilter'));
            M.grob.parset = uicontrol('style','edit','max',2, ...
                'backgroundcolor','white','horizontalalignment','left', ...
                'callback',@(u,evt)parest_manage(M,'changepar'));
            M.grob.paractive = uicontrol('style','checkbox','string','active', ...
                'backgroundcolor',bgcolor, ...
                'callback',@(u,evt)parest_manage(M,'active'));
            % (estimate)
            M.grob.huestimate = uicontrol('string','ESTIMATE', ...
                'callback',@(u,evt)res_estimate(M,'flow'));
            M.grob.huestimateall = uicontrol('string','EST. ALL', ...
                'callback',@(u,evt)res_estimate(M,'flowall'));

            % SECTIONS
            % (manage parameters list)
            M.grob.sectionlist = uicontrol('style','listbox', ...
                'callback',@(u,evt)sections_manage(M,'select',get(u,'value')));
            M.grob.sectionlisttop  = uicontrol('string','move to top', ...
                'callback',@(u,evt)sections_manage(M,'top'));
            M.grob.sectionlistrm   = uicontrol('string',' remove ', ... MATLAB BUG! 'remove' produces blank display! use ' remove ' instead
                'callback',@(u,evt)sections_manage(M,'remove'));
            M.grob.sectionlistsave = uicontrol('string','add', ...
                'callback',@(u,evt)sections_manage(M,'add'));
            M.grob.sectionlistreplace = uicontrol('string','replace', ...
                'callback',@(u,evt)sections_manage(M,'replace'));
            % (define parameters)
            M.grob.sectionid = uicontrol('style','text');
            M.grob.sectionpar = uicontrol('style','edit','max',2, ...
                'backgroundcolor','white','horizontalalignment','left', ...
                'callback',@(u,evt)sections_manage(M,'changesectionpar'));
            M.grob.sectionestpar = uicontrol('style','edit','max',2, ...
                'backgroundcolor','white','horizontalalignment','left', ...
                'callback',@(u,evt)sections_manage(M,'changeestpar'));
            M.grob.sectionactive = uicontrol('style','checkbox','string','active', ...
                'backgroundcolor',bgcolor, ...
                'callback',@(u,evt)sections_manage(M,'active'));
            % (estimate)
            M.grob.sectionestimate = uicontrol('string','ESTIMATE', ...
                'callback',@(u,evt)res_estimate(M,'section'));
            M.grob.sectionestimateall = uicontrol('string','EST. ALL', ...
                'callback',@(u,evt)res_estimate(M,'sectionall'));
           
            % TRIALS
            M.grob.hatrials = axes('fontsize',10);
            M.misc.hutrial = fn_stepper('min',0,'max',0,'step',1, ...
                'format','%i','coerce',true,'value',0, ...
                'callback',@(u,evt)trials_setcurrent(M,get(u,'value')));
            M.grob.hutrial = M.misc.hutrial.hpanel; % need the grob value to be a Matlab graphic object, not a fn_stepper object
            
            % RESULTS
            for i=1:7, M.grob.hares(i) = axes('units','pixel','fontsize',10); end
            M.grob.hulinshift = uicontrol('style','popupmenu', ...
                'string',{'no shift', ...
                '1pix','2pix','3pix','4pix','5pix','6pix','7pix','8pix','9pix','10pix', ...
                '-1pix','-2pix','-3pix','-4pix','-5pix','-6pix','-7pix','-8pix','-9pix','-10pix'}, ...
                'value',1, ...
                'callback',@(u,evt)res_display(M,[],true));
        end
        function init_menus(M)
            F = fieldnames(M.menus);
            for i=1:length(F)
                delete(M.menus.(F{i}))
            end
            init_menus@interface(M)
            m = M.menus.interface;
            uimenu(m,'label','re-display','separator','on', ...
                'callback',@(u,evt)fastflow_display(M))
            uimenu(m,'label','object in basespace', ...
                'callback',@(u,evt)assignin('base','M',M))
            uimenu(m,'label','new wide field experiment','separator','on', ...
                'callback',@(u,evt)fastflow_manage(M,'newwidefield'))
            uimenu(m,'label','new line scan experiment', ...
                'callback',@(u,evt)fastflow_manage(M,'newlinescan'))
            uimenu(m,'label','new fake', ...
                'callback',@(u,evt)fastflow_manage(M,'newfake'))
            uimenu(m,'label','open','separator','on', ...
                'callback',@(u,evt)fastflow_manage(M,'open'))
            uimenu(m,'label','save', ...
                'callback',@(u,evt)fastflow_manage(M,'save'))
            uimenu(m,'label','registration & average movie settings','separator','on', ...
                'callback',@(u,evt)fastflow_manage(M,'settings'))
            uimenu(m,'label','interp+vol','separator','on', ...
                'callback',@(u,evt)fastflow_manage(M,'interp+vol'))
            uimenu(m,'label','analysis', ...
                'callback',@(u,evt)fastflow_manage(M,'analysis'))
            uimenu(m,'label','results vol','separator','on', ...
                'callback',@(u,evt)fastflow_manage(M,'resultsvol'))
            uimenu(m,'label','vessel color', ...
                'callback',@(u,evt)fastflow_manage(M,'vesselcolor'))
            
            % vessel segmentation
            vessels_menu(M)            
            
            % estimation parameters
            parest_menu(M)
            
            % sections parameters
            sections_menu(M)
            
            % result display
            res_menu(M)
        end
        function misc_init(M)            
            % figure
            set(M.hf,'windowbuttondownfcn',@(hf,evnt)misc_axesclick(M,gca))
            % structure 'haschanged'
            M.misc.haschanged = struct('ktrial',false,'trials',false, ...
                'kedge',false,'algo',false,'algos',false,'fastflow',false, ...
                'section',false,'sections',false);
            % action button
            hua = M.grob.actionbutton;
            m = uicontextmenu('parent',M.hf);
            set(hua,'enable','off', ...
                'uicontextmenu',m,'callback','fast_actionbutton')
            m1=uimenu(m,'label','enable action button','checked',get(hua,'enable'), ...
                'callback',@enabledisable);
            uimenu(m,'label','edit action code', ...
                'callback',@(hu,e)edit('fast_actionbutton'))
            function enabledisable(hu,e) %#ok<INUSD>
                val = get(u,'enable');
                val = fn_switch(val,'switch');
                set(u,'enable',val)
                set(m1,'checked',val)
            end
            % result display holding
            M.misc.holdresdisplay = false;
        end
    end
    
    % Content
    methods
        function fastflow_manage(M,flag)
            switch flag
                case 'start'
                    % check if we display an image
                    M.im.doim = status(M.S,'movie');
                    % display the fastflow object
                    fastflow_display(M)
                case {'newwidefield' 'newlinescan' 'newfake'}
                    fastflow_cleanup(M), drawnow
                    closeproject(M.S)
                    M.S = fastflow;
                    switch flag
                        case 'newwidefield'
                            inputwidefield(M.S)
                        case 'newlinescan'
                            inputlinescan(M.S)
                        case 'newfake'
                            inputfake(M.S)
                    end
                    % check if we display an image
                    M.im.doim = status(M.S,'movie');
                    if M.im.doim
                        M.S.segmentpar = M.options.segment; % default segmentation parameters
                    end
                    % display the fastflow object
                    fastflow_display(M)
                case 'open'
                    % close some figures related to the current fastflow
                    % object
                    closeproject(M.S)
                    % file name
                    fname = fn_getfile('*_batch.mat', ...
                        'Select existing batch data file');
                    if ~fname, return, end
                    % load the fastflow object
                    M.S = fn_loadvar(fname);
                    % update immediately savedir and nickname based on file
                    % name
                    [rep nick] = fileparts(fname);
                    M.S.savedir = [rep '/'];
                    M.S.nickname = strrep(strrep(nick,'_batch',''),'_light','');
                    % handle old versions
                    if isempty(M.S.segmentpar), M.S.segmentpar = M.options.segment; end
                    % check data and resample directories
                    updatefolders(M.S)
                    % check if we display an image
                    M.im.doim = status(M.S,'movie');
                    % display the fastflow object
                    fastflow_display(M)
                case 'save'
                    deletefiles = cleanfiles(M.S);
                    if ~isempty(deletefiles) 
                        if ~fn_reallydlg( ...
                                'By saving, the following files will be deleted:', ...
                                deletefiles, 'Do you want to proceed?');
                            return
                        end
                        deletefiles = strcat(M.S.savedir,deletefiles);
                        delete(deletefiles{:});
                    end
                    saveproject(M.S)
                case 'settings'
                    setparameters(M.S)
                case 'interp+vol'
                    % save first
                    fastflow_manage(M,'save')
                    analysis(M.S,'interp+vol')
                case 'analysis'
                    % save first
                    fastflow_manage(M,'save')
                    analysis(M.S)
                case 'resultsvol'
                    resultsvol(M.S)
                case 'vesselcolor'
                    fast_color(M.S)
                otherwise
                    error('unknown flag ''%s''',flag)
            end
        end
        function fastflow_cleanup(M)
            set(M.misc.hutrial,'min',0,'max',0,'value',0)
            M.ktrial = [];
            set(M.hf,'name','FASTFLOW')
            vessels_cleanup(M)
            trials_cleanup(M)
        end
        function fastflow_display(M)
            % clean-up first
            fastflow_cleanup(M)
            % show/hide specific menus and controls
            if M.im.doim
                onoffim = 'on';
                onoffscan = 'off';
            else
                onoffim = 'off';
                onoffscan = 'on';
            end
            set([M.vsl.menuitems M.im.objects],'visible',onoffim)
            set([M.scan.menuitems M.scan.objects],'visible',onoffscan)
            % exit if there is no fastflow content
            if ~status(M.S,'base')
                return
            end
            % modify the resdisp structures to take into account change in
            % S.parameters.fastmain.resultavgtrials and S.cond
            for i=1:length(M.res)
                M.res(i).s = fast_resdispselect(M.S,M.res(i).s,'correction');
            end
            % figure name
            set(M.hf,'name',['FASTFLOW - ' M.S.nickname])
            % zooming enabled
            fn_imvalue %image
            % current trial
            if M.im.doim
                set(M.misc.hutrial,'min',1,'max',M.S.nexp)
                trials_setcurrent(M,1)
                % display trials
                trials_display(M)
            else
                set(M.misc.hutrial,'min',1,'max',Inf)
                trials_setcurrent(M,1)
                trials_cleanup(M)
            end
            % display vessels
            M.misc.haschanged.fastflow = true;
            vessels_display(M)
            % display the color table for groups
            parest_displaygroups(M)
            % current vessel
            if ~isempty(M.S.edges)
                % this will display vessel, estimation parameters and
                % results
                vessels_setcurrent(M,1)
            else
                % update result display (mostly erase graphs!)
                M.misc.haschanged.kedge = true;
                M.misc.haschanged.algo = true;
                res_display(M)
            end
        end
    end
    
    % Vessels
    methods
        % first initialization
        function vessels_init(M)
            % colormap gray for the whole figure
            colormap(M.grob.hatrials,'gray')
            % scan
            M.scan.order = [];
            M.scan.kedges = [];
            M.scan.menu = uicontextmenu('parent',M.hf);
            m = M.scan.menu;
            uimenu(m,'label','set active/inactive', ...
                'callback',@(hu,evnt)scan_manage(M,'toggleactive'))
            uimenu(m,'label',' remove ', ... % MATLAB BUG!!
                'callback',@(hu,evnt)scan_manage(M,'remove'))
            set(M.grob.scanlist,'UIContextMenu',m)
            % image
            M.im.showim = false;
            M.im.ylimlisten = [];
            M.im.timer = [];
            % vessels. Be careful: M.vsl.showwidth, showwidthmenu,
            % zoommode, zoommodemenu already initialized in vessels_menu 
            M.vsl.menu = uicontextmenu('parent',M.hf);
            m = M.vsl.menu;
            uimenu(m,'label','set active/inactive', ...
                'callback',@(hu,evnt)vessels_manage(M,gco,'toggleactive'))
            uimenu(m,'label','upside-down','separator','on', ...
                'callback',@(hu,evnt)vessels_manage(M,gco,'upsidedown'))
            uimenu(m,'label','re-segment edge', ...
                'callback',@(hu,evnt)vessels_manage(M,gco,'resegment'))
            uimenu(m,'label','delete edge','separator','on', ...
                'callback',@(hu,evnt)vessels_manage(M,gco,'delete'))
        end
        function vessels_menu(M)
            if ~isfield(M.vsl,'showwidth'), M.vsl.showwidth=true; end
            if ~isfield(M.vsl,'zoommode'), M.vsl.zoommode=false; end
            if ~isfield(M.vsl,'dotrialframe'), M.im.dotrialframe=false; end
            if ~isfield(M.im,'doim'), M.im.doim = false; end
            mv = [];
            mf = [];
            m = uimenu(M.hf,'label','Vessels');
            mv(end+1) = uimenu(m,'label','movie +', ...
                'callback',@(u,evt)im_movie(M));
            mv(end+1) = uimenu(m,'label','show width','separator','on', ...
                'checked',fn_switch(M.vsl.showwidth), ...
                'callback',@(u,evt)vessels_manage(M,[],'showwidth',~M.vsl.showwidth));
            M.vsl.showwidthmenu = mv(end); 
            mv(end+1) = uimenu(m,'label','zoom mode', ...
                'checked',fn_switch(M.vsl.zoommode), ...
                'callback',@(u,evt)vessels_manage(M,[],'zoommode',~M.vsl.zoommode));
            M.vsl.zoommodemenu = mv(end); 
            mv(end+1) = uimenu(m,'label','first frame of current trial', ...
                'checked',fn_switch(M.im.dotrialframe), ...
                'callback',@(u,evt)vessels_manage(M,[],'dotrialframe',~M.im.dotrialframe));
            M.im.dotrialframemenu = mv(end);
            uimenu(m,'label','all active','separator','on', ...
                'callback',@(u,evt)vessels_manage(M,[],'allactive',true))
            uimenu(m,'label','all inactive', ...
                'callback',@(u,evt)vessels_manage(M,[],'allactive',false))
            mv(end+1) = uimenu(m,'label','edit segmentation parameters...','separator','on', ...
                'callback',@(u,evt)vessels_segmentpar(M,'edit'));
            mv(end+1) = uimenu(m,'label','go back to default', ...
                'callback',@(u,evt)vessels_segmentpar(M,'getdefault'));
            mv(end+1) = uimenu(m,'label','make current default', ...
                'callback',@(u,evt)vessels_segmentpar(M,'setdefault'));
            mv(end+1) = uimenu(m,'label','delete current','separator','on', ...
                'callback',@(u,evt)vessels_manage(M,'current','delete'));
            M.menus.vessels = m;
            M.vsl.menuitems = mv;
            M.scan.menuitems = mf;
            if M.im.doim, set(mf,'visible','off'), else set(mv,'visible','off'), end
        end
        % init/endup display
        function vessels_display(M)
            ha = M.grob.haim; cla(ha)
            M.vsl.hl = []; 
            ne = length(M.S.edges);
            % list of vessels ('scans') or image?
            if M.im.doim
                % image
                M.im.him = imagesc(M.S.CS','parent',ha);
                set(ha,'xtick',[],'ytick',[])
                axis(ha,'image')
                im_toggle(M)
                % vessels
                M.vsl.hl = zeros(ne,5);
                for i=1:ne
                    vsl_display(M,i) % draw vessel i and store handles in M.vsl.hl(i,:)
                end
            else
                % display list of vessels
                ne = length(M.S.edges);
                list = cell(1,ne);
                for i=1:ne
                    e = M.S.edges(i);
                    switch e.interpflag
                        case 'movie'
                            str = e.filebase;
                        case 'scan'
                            str = e.filepar.filebase;
                        case 'fake'
                            s = e.fakepar;
                            if s.step>=10,      step = [num2str(floor(s.step))      '.']; else step = num2str(s.step,     '%.1f'); end
                            if s.shotnoise>=10, noise= [num2str(floor(s.shotnoise)) '.']; else noise= num2str(s.shotnoise,'%.1f'); end
                            if mod(s.density,.1), dens = ['den .' num2str(floor(s.density*100),'%.2i')]; else dens = ['dens .' num2str(s.density*10)]; end
                            str = sprintf( ...
                                ['step %s - %s - diam %.2i - jit %.2i - ' ...
                                'noise %s [%.3ix%.3ix%.3i]'], ...
                                step,dens,floor(s.diam),floor(s.jitter), ...
                                noise,e.np,e.nt,e.nexp);
                    end
                    brackets=fn_switch(e.active,'[]','()');
                    numstr = sprintf('%c%i%c',brackets(1),e.flagnumber,brackets(2));
                    list{i} = [numstr ' ' str];
                end
                set(M.grob.scanlist,'string',list,'value',M.scan.kedges)
            end
            % button for vessel number fast selection
            set(M.misc.huvessel,'min',1,'max',ne)
        end
        function vessels_cleanup(M)
            % unselect current vessel
            vessels_unselect(M)
            set(M.misc.huvessel,'min',0,'max',0)
            % clean-up image display
            cla(M.grob.haim)
            set(M.grob.scanlist,'string','','value',1)
            M.vsl.hl = [];
            delete(M.im.ylimlisten)
            M.im.ylimlisten = [];
        end
        % set/unset current
        function vessels_unselect(M)
            if isempty(M.kedge), return, end
            if M.im.doim
                % make line normal (thinner)
                try set(M.vsl.hl(M.kedge,:),'linewidth',1), end
            end
            % no current edge
            M.kedge = [];
            M.e = fast_edge.empty;
            M.misc.haschanged.kedge = true;
            set(M.misc.huvessel,'value',0)
            % no pixel shift in image display
            set(M.grob.hulinshift,'value',1)
            % update parameters and result display
            parest_cleanup(M)
            sections_cleanup(M)
            res_display(M)
        end
        function vessels_setcurrent(M,k)
            if k==0, return, end 
            % unmark current first (this is only part of 'vessels_cleanup'
            % code)
            if M.im.doim
                try set(M.vsl.hl(M.kedge,:),'linewidth',1), end
            end
            set(M.grob.hulinshift,'value',1)
            M.kflow = [];
            % properties
            M.kedge = k;
            M.e = M.S.edges(M.kedge);
            M.misc.haschanged.kedge = true;
            set(M.misc.huvessel,'value',M.kedge)
            % display
            if M.im.doim
                % make line thicker
                set(M.vsl.hl(M.kedge,:),'linewidth',3)
            else
                % select edge (and allow only one selection)
                set(M.grob.scanlist,'value',k)
            end
            if M.vsl.zoommode, vsl_zoom(M,k), end
            drawnow
            % display vessel estimations and update result display
            M.misc.holdresdisplay = true; % parest AND sections should be updated before display update
            parest_display(M)
            if M.im.doim, sections_display(M), end
            M.misc.holdresdisplay = false;
            res_display(M)
        end
        % image functions
        function im_toggle(M)
            global Y
            if ~status(M.S,'base'), return, end
            ha  = M.grob.haim;
            % clean
            delete(M.im.ylimlisten)
            M.im.ylimlisten = [];
            if ~isempty(M.im.timer)
                stop(M.im.timer)
                delete(M.im.timer)
                M.im.timer = [];
            end
            kvisible = fn_switch(M.vsl.showwidth,1:5,[1 5]);
            if ~isempty(M.vsl.hl), set(M.vsl.hl(:,kvisible),'visible','on'), end
            % choose what to show
            label = get(get(M.grob.im_toggle,'selectedobject'),'string');
            switch label
                case 'vessels'
                    im_displayframe(M,true)
                case 'energy'
                    maskcol = M.S.mask' / (2*std(M.S.mask(:)));
                    maskcol = 1+max(0,min(255,round((maskcol+1)/2*255)));
                    cm = signcheck(256);
                    [ni nj] = size(maskcol);
                    maskcol = reshape(cm(maskcol(:),:),[ni nj 3]);
                    set(M.im.him,'cdata',maskcol);
                case 'movie'
                    % hide vessels
                    set(M.vsl.hl,'visible','off')
                    % compute frame average
                    loadtrial(M.S,M.ktrial)
                    M.im.fravg = mean(Y,3);
                    if ~isa(Y,'double'), M.im.fravg = single(M.im.fravg); end
                    % cliping 
                    set(ha,'clim',[-1 1]*.02)
                    % start movie
                    M.im.kframe = 1;
                    M.im.timer = timer( ...
                        'ExecutionMode',    'fixedSpacing', ...
                        'Period',           .066, ...
                        'TimerFcn',         @(u,e)im_movieframe());
                    start(M.im.timer)
            end
            function im_movieframe
                frame = fn_float(Y(:,:,M.im.kframe)) ./ M.im.fravg;
                frcut = frame(6:end-5,6:end-5);
                m = mean(frcut(~isnan(frcut)));
                frame(isnan(m)) = m;
                frame = frame - m;
                set(M.im.him,'cdata',frame')
                M.im.kframe = 1+mod(M.im.kframe,M.S.nt-1);
            end
        end
        function im_movie(M)
            global Y
            if ~M.im.doim, error programming, end
            loadtrial(M.S,M.ktrial)
            fn_movie(Y)
        end
        function im_displayframe(M,doautoclip)
            global Y
            if nargin<2, doautoclip=false; end
            if M.im.dotrialframe
                doloadtrial = strcmp(M.S.parameters.registration,'none') ...
                    || isequal(M.S.ycontains,[M.ktrial 1]);
                if ~doloadtrial
                    try
                        % read only the first frame
                        [dum b] = fileparts(M.S.files(M.ktrial,:)); %#ok<ASGLU>
                        y = fast_load2bytes([M.S.resampledir b '.ff'],1);
                    catch
                        doloadtrial = true;
                    end
                end
                if doloadtrial
                    loadtrial(M.S,1)
                    y = Y(:,:,1);
                end
            else
                y = M.S.CS;
            end
            set(M.im.him,'cdata',y');
            ha = M.grob.haim;
            M.im.ylimlisten = addlistener(ha,'YLim', ...
                'PostSet',@(m,evnt)im_autoclip());
            if doautoclip, im_autoclip, end
            % special: adjust clipping after zooming
            function im_autoclip
                ax = axis(ha);
                try %#ok<TRYNC>
                    ii = round(ax(1)):round(ax(2)-1e-6);
                    jj = round(ax(3)):round(ax(4)-1e-6);
                    x = M.S.CS(ii,jj);
                    clip = [min(x(:)) max(x(:))];
                    if diff(clip)==0, clip = clip+[-1 1]; end
                    set(ha,'clim',clip)
                end
            end
        end
        % segmentation
        function vessels_selectpoly(M,ha,flag)
            p = get(ha,'currentpoint'); p=p(1,1:2)';
            % input flag
            if nargin<3
                switch get(M.hf,'selectiontype')
                    case 'open'             % start/finish new poly selection
                        if isempty(M.polysel)
                            flag = 'new';
                        else
                            flag = 'end';
                        end
                    case 'normal'           % add a point to the poly
                        if isempty(M.polysel), return, end
                        flag = 'point';
                    case {'alt' 'extend'}   % cancel the current selection
                        if isempty(M.polysel), return, end
                        flag = 'cancel';
                end
            end
            % action
            switch flag
                case 'new'
                    M.polysel = struct;
                    M.polysel.poly = [p p];
                    M.polysel.hl = line(p(1),p(2),'parent',ha, ...
                        'color',[.5 .5 .5]);
                    set(M.hf,'windowbuttonmotionfcn',@(h,evnt)vessels_selectpoly(M,ha,'move'))
                case 'move'
                    M.polysel.poly(:,end) = p;
                    set(M.polysel.hl,'xdata',M.polysel.poly(1,:),'ydata',M.polysel.poly(2,:))
                case 'point'
                    M.polysel.poly(:,end+1) = p;
                    set(M.polysel.hl,'xdata',M.polysel.poly(1,:),'ydata',M.polysel.poly(2,:))
                case 'cancel'
                    delete(M.polysel.hl)
                    M.polysel = [];
                    set(M.hf,'windowbuttonmotionfcn','')
                case 'end'
                    poly = M.polysel.poly(:,1:end-1);
                    if size(poly,2)==1, return, end % only one point selected
                    delete(M.polysel.hl)
                    M.polysel = [];
                    set(M.hf,'windowbuttonmotionfcn','')
                    % automatic segmentation
                    s = fast_segment(M.S.mask,poly,M.S.segmentpar,ha);
                    % store new vessel
                    if ~status(M.S,'movie'), error programming, end
                    edge = fast_edge(M.S,s.snake,s.halfd,s.dx);
                    % set the estimation parameters according to
                    % default
                    parest_groupdefault(M,'setedge',edge)
                    k = length(M.S.edges)+1;
                    set(M.misc.huvessel,'max',k) % update the display for selection of vessel number
                    M.S.edges(:,k) = edge;
                    vsl_display(M,k)
                    vessels_setcurrent(M,k)
            end
        end
        function vessels_segmentpar(M,flag)
            switch flag
                case 'edit'
                    p = M.S.segmentpar;
                    a = inputdlg('segmentation parameters','edit parameters', ...
                        length(p)+2,{strvcat(p{:})}); %#ok<VCAT>
                    if ~isempty(a)
                        M.S.segmentpar = cellstr(a{1});
                    end
                case 'getdefault'
                    M.S.segmentpar = M.options.segment;
                case 'setdefault'
                    M.options.segment = M.S.segmentpar;
                    saveoptions(M)
            end
        end
        % vessel functions
        function vsl_display(M,i)
            ha = M.grob.haim;
            ei = M.S.edges(1,i);
            p = ei.snake(:,ceil(length(ei.snake)/3));
            hl = zeros(1,5);
            hl(1)=line(ei.snake(1,:),ei.snake(2,:), ...
                'parent',ha,'hittest','off', ...
                'userdata',i, ...
                'buttondownfcn',@(hp,evnt)vsl_click(M,get(hp,'userdata')), ...
                'uicontextmenu',M.vsl.menu);
            hl(2)=line(ei.tubea(1,:),ei.tubea(2,:),'linestyle','--', ...
                'parent',ha,'hittest','off');
            hl(3)=line(ei.tubeb(1,:),ei.tubeb(2,:),'linestyle','--', ...
                'parent',ha,'hittest','off');
            hl(4)=line(p(1),p(2),'marker','o','markersize',10, ...
                'parent',ha,'hittest','on', ...
                'userdata',i, ...
                'buttondownfcn',@(hp,evnt)vsl_click(M,get(hp,'userdata')), ...
                'uicontextmenu',M.vsl.menu);
            hl(5)=text(double(ei.snake(1,1)),double(ei.snake(2,1)),num2str(i), ...
                'parent',ha,'hittest','off');
            if ~M.vsl.showwidth
                set(hl(2:4),'visible','off')
                set(hl(1),'hittest','on')
            end
            M.vsl.hl(i,:) = hl;
            vsl_color(M,i)
        end
        function vsl_color(M,i)
            [col col0] = getcolor(M.S.edges(1,i));
            hl = M.vsl.hl(i,:);
            set(hl(4),'col',col0)         % color according to flag for the circle
            set(hl([1:3 5]),'col',col)    % same for the rest, but black if edge not active
        end
        function vsl_click(M,i)
            if ~strcmp(get(M.hf,'selectiontype'),'normal'), return, end
            vessels_setcurrent(M,i);
        end
        function vsl_zoom(M,i)
            if nargin<2, i=M.kedge; end
            if isempty(i), return, end
            ei = M.S.edges(i);
            p = mean(ei.snake,2);
            nside = max(max(ei.snake(:,[1 end])-[p p]));
            nside = ceil(nside*[3 2.3]);
            p = max([1+nside(1);1+nside(2)],min([M.S.nx-nside(1);M.S.ny-nside(2)],round(p)));
            newax = [p(1)+[-1 1]*nside(1) p(2)+[-1 1]*nside(2)];
            fn_imvalue('chgxy',newax,M.grob.haim)
        end
        function vsl_clear(M,ivessel)
            % files to delete
            ei = M.S.edges(1,ivessel);
            if ~isempty(ei.flow) && ~isempty([ei.flow.files]) ...
                    && ~fn_reallydlg({'Estimation results might exist,' ...
                    'are you sure you want to loose track of them?'});
                return
            end
            % remove line display
            delete(M.vsl.hl(ivessel,:))
            % do not remove the fast_edge element from S yet, because we
            % don't know if it has to be a deletion or a replacement, BUT
            % THE CALLING FUNCTION MUST DO IT
        end
        % vessel menu actions
        function vessels_manage(M,hp,flag,value)
            if ~isempty(hp)
                if ischar(hp)
                    ivessel = M.kedge;
                elseif ishandle(hp)
                    ivessel = get(hp,'userdata');
                else
                    error('argument')
                end
                ei = M.S.edges(ivessel);
            end
            switch flag
                % CURRENT VESSEL
                case 'toggleactive'
                    [ei.active] = deal(~ei.active);
                    vsl_color(M,ivessel)
                case 'upsidedown'
                    % delete files, erase vessel display
                    vsl_clear(M,ivessel)
                    % create a new vessel with points inverted
                    edge = fast_edge(M.S,ei.snake(:,end:-1:1), ...
                        ei.halfd(:,end:-1:1),ei.dx);
                    % store new vessel
                    if ~status(M.S,'movie'), error programming, end
                    M.S.edges(:,ivessel) = edge;
                    parest_groupdefault(M,'setedge',edge)
                    % re-display 
                    vsl_display(M,ivessel)
                    vessels_setcurrent(M,ivessel)
                case 'resegment'
                    % delete files, erase vessel display
                    vsl_clear(M,ivessel)
                    % automatic segmentation
                    s = fast_segment(M.S.mask,ei.snake, ...
                        M.S.segmentpar,M.grob.haim);
                    % store new vessel
                    if ~status(M.S,'movie'), error programming, end
                    edge = fast_edge(M.S,s.snake,s.halfd,s.dx);
                    M.S.edges(:,ivessel) = edge;
                    parest_groupdefault(M,'setedge',edge)
                    % re-display 
                    vsl_display(M,ivessel)
                    vessels_setcurrent(M,ivessel)
                case 'delete'
                    % delete files, erase vessel display
                    if ~isempty(ei.flow) && ~isempty([ei.flow.files]) ...
                            && ~fn_reallydlg({'Estimation results might exist,' ...
                            'are you sure you want to loose track of them?'});
                        return
                    end
                    % remove edge
                    if M.kedge==ivessel, vessels_unselect(M), end
                    delete(M.vsl.hl(ivessel,:))
                    M.vsl.hl(ivessel,:) = [];            
                    M.S.edges(:,ivessel) = [];
                    % need to update some numbers!!!
                    ne = length(M.S.edges);
                    for i=1:ne
                        set(M.vsl.hl(i,[1 4]),'userdata',i)
                        set(M.vsl.hl(i,5),'string',num2str(i))
                    end
                    set(M.misc.huvessel,'min',1,'max',ne)
                    if M.kedge>ivessel
                        M.kedge=M.kedge-1; 
                        set(M.misc.huvessel,'value',M.kedge)
                    end 
                % ALL VESSELS
                case 'showwidth'
                    M.vsl.showwidth = value;
                    set(M.vsl.showwidthmenu,'checked',fn_switch(value))
                    if value
                        set(M.vsl.hl(:,2:4),'visible','on')
                        set(M.vsl.hl(:,1),'hittest','off')
                    else
                        set(M.vsl.hl(:,2:4),'visible','off')
                        set(M.vsl.hl(:,1),'hittest','on')
                    end
                case 'allactive'
                    [M.S.edges.active] = deal(value);
                    if M.im.doim
                        for i=1:length(M.S.edges), vsl_color(M,i), end
                    else
                        scan_manage(M,'color')
                    end
                case 'zoommode'
                    M.vsl.zoommode = value;
                    set(M.vsl.zoommodemenu,'checked',fn_switch(value))
                    if value, vsl_zoom(M), end
                % IMAGE
                case 'dotrialframe'
                    M.im.dotrialframe = value;
                    set(M.im.dotrialframemenu,'checked',fn_switch(value))
                    
            end
        end
        % list of edges
        function scan_manage(M,flag,value)
            switch flag
                case 'color'
                    vessels_display(M)
                case 'newfake'
                    if isempty(M.kedge)
                        par0=fast_fakedata('par'); 
                    else
                        par0=M.e.fakepar; 
                    end
                    par = fast_fakedata('userpar',par0);
                    if isempty(par), return, end
                    ke = length(M.S.edges)+1;
                    M.S.edges(ke) = fast_edge(M.S,'fake',par);
                    vessels_display(M)
                    vessels_setcurrent(M,ke)
                case 'remove'
                    value = M.scan.kedges;
                    if isempty(value), return, end
                    flows = [M.S.edges(value).flow];
                    if ~isempty([flows.files]) && ~fn_reallydlg( ...
                            'Estimation results might exist,', ...
                            'are you sure you want to loose track of them?')
                        return
                    end
                    M.S.edges(value) = [];
                    vessels_display(M)
                    ne = length(M.S.edges);
                    if ne
                        M.scan.kedges = min(value(1),ne);
                        vessels_setcurrent(M,M.scan.kedges)
                    else
                        set(M.grob.scanlist,'value',[])
                        M.scan.kedges = [];
                        vessels_unselect(M)
                    end
                case 'reorder'
                    % color the button to remember that we are in 'reorder'
                    % mode when selecting in the list of edges
                    if isempty(M.scan.order)
                        set(M.grob.scanorder,'backgroundcolor','b')
                        set(M.grob.scanlist,'value',[])
                        M.scan.kedges = [];
                        vessels_unselect(M)
                        M.scan.order = 0; % number of parameters ordered
                    else
                        % finish ordering
                        ne = length(M.S.edges);
                        set(M.grob.scanlist,'value',1:ne), pause(.2)
                        set(M.grob.scanlist,'value',[])
                        set(M.grob.scanorder,'backgroundcolor','default')
                        M.scan.order = [];
                    end
                case 'select'
                    if ~isempty(M.scan.order)
                        % re-ordering mode
                        ord = M.scan.order;
                        ne = length(M.S.edges);
                        if isscalar(value) && value>ord
                            perm = [1:ord value setdiff(ord+1:ne,value)];
                            ord = ord+1;
                        elseif length(value)>ord && all(value==1:length(value))
                            perm = 1:ne;
                            ord = length(value);
                        else
                            % bad selection
                            set(M.grob.scanlist,'value',1:ord)
                            return
                        end
                        list = get(M.grob.scanlist,'string');
                        list = list(perm);
                        M.S.edges = M.S.edges(perm);
                        set(M.grob.scanlist,'string',list,'value',1:ord)
                        if ord>=ne-1
                            % with this last element, ordering is finished
                        	set(M.grob.scanlist,'value',1:ne-1), pause(.2)
                        	set(M.grob.scanlist,'value',1:ne),   pause(.2)
                        	set(M.grob.scanlist,'value',[])
                            set(M.grob.scanorder,'backgroundColor','default')
                            M.scan.order = [];
                        else
                        	set(M.grob.scanlist,'value',1:ord)
                            M.scan.order = ord;
                        end
                    elseif isscalar(value)
                        % select one vessel
                        M.scan.kedges = value;
                        vessels_setcurrent(M,value)
                    else
                        % select several vessels
                        M.scan.kedges = value;
                        vessels_unselect(M)
                    end
                case 'toggleactive'
                    if isempty(M.scan.kedges), return, end
                    edges = M.S.edges(M.scan.kedges);
                    b = ~all([edges.active]);
                    [edges.active] = deal(b);
                    vessels_display(M)
                otherwise
                    error('unknown flag ''%s''',flag)
            end
        end
    end
        
    % Estimation parameters
    methods
        % first initialization
        function parest_init(M)
            M.parest = struct;
            M.parest.methodlist = {'track' 'radon' 'gabor'};
            set(M.grob.parmethod,'string',M.parest.methodlist, ...
                'callback',@(u,evt)parest_manage(M,'method',M.parest.methodlist(get(u,'value'))))
            M.parest.controls = [M.grob.info(2) M.grob.haflags M.grob.parlist ...
                M.grob.parlisttop M.grob.parlistrm M.grob.parlistsave M.grob.parlistreplace ...
                M.grob.parmethod M.grob.parid M.grob.parfilter M.grob.parset ...
                M.grob.paractive M.grob.huestimate M.grob.huestimateall];             
        end
        function parest_menu(M)
            m = uimenu(M.hf,'label','Parameters');
            uimenu(m,'label','reorder [track-radon-gabor]', ...
                'callback',@(u,evt)parest_manage(M,'all_reorder_trg'))
            uimenu(m,'label','reorder [radon-track-gabor]', ...
                'callback',@(u,evt)parest_manage(M,'all_reorder_rtg'))
            uimenu(m,'label','add current to all vessels','separator','on', ...
                'callback',@(u,evt)parest_manage(M,'all_addcurrent'))
            uimenu(m,'label','change for all vessels...', ...
                'callback',@(u,evt)parest_changeallvessels(M))
            uimenu(m,'label','make currents group default','separator','on', ...
                'callback',@(u,evt)parest_groupdefault(M,'setdefault'))
            uimenu(m,'label','... and add to group members', ...
                'callback',@(u,evt)parest_manage(M,'group_add'))
            uimenu(m,'label','... and replace group members', ...
                'callback',@(u,evt)parest_manage(M,'group_replace'))
            uimenu(m,'label','edit groups...', ...
                'callback',@(u,evt)parest_editgroups(M))
            M.menus.parest = m;
        end
        % display
        function parest_cleanup(M)
            % no selected algorithm
            M.kflow = [];
            M.f = fast_flow.empty;
            % empty algo list
            set(M.grob.parlist,'string','','value',0)
            % unmark group
            set(M.parest.hlgroupmark,'visible','off')
        end
        function parest_display(M)
            if isempty(M.kedge), error programming, end
            nf = length(M.e.flow);
            % mark group to which edges belong
            set(M.parest.hlgroupmark,'visible','on','xdata',getflagnumber(M.e))
            % any existing flow? if not, create one
            if nf>0
                M.parest.unsaved = false;
                if isempty(M.kflow)
                    M.kflow = 1;
                else
                    M.kflow = min(M.kflow,nf);
                end
                M.f = M.e.flow(M.kflow); 
            else
                M.parest.unsaved = true;
                method = M.parest.methodlist{get(M.grob.parmethod,'value')};
                estpar = defaultpar(M,method);
                M.f = fast_flow(M.e,method,M.options.filterpar,estpar);
            end
            drawnow
            % then, display everything
            parest_displaypar(M)
        end
        function parest_displaypar(M)
            % list of algos present in current edge
            if isempty(M.e.flow)
                % MATLAB BUG!!!
                algolist = {};
            else
                algolist = {M.e.flow.algo};
            end
            %disp(algolist)
            if M.parest.unsaved
                algolist{end+1} = '(not saved)';
                klist = length(algolist);
            else
                klist = M.kflow;
            end
            set(M.grob.parlist,'string',algolist,'value',klist)
            % method
            kmethod = find(strcmp(M.parest.methodlist,M.f.method));
            set(M.grob.parmethod,'value',kmethod)
            % ID
            set(M.grob.parid,'string',M.f.id)
            % filtering parameters
            set(M.grob.parfilter,'string',M.f.filterpar)
            % estimation parameters
            set(M.grob.parset,'string',fn_struct2str(M.f.estpar))
            % active flag
            set(M.grob.paractive,'value',M.f.active)
            % update result display
            res_display(M)
        end
        function parest_displaygroups(M)
            flaglist = getflaglist(M.S);
            cols = shiftdim(cat(1,flaglist.color),-1);
            set(M.hf,'tag','haflags') % prevents fn_imvalue to have effect on this graph
            image(cols,'parent',M.grob.haflags, ...
                'hittest','on','buttondownfcn',@(u,e)parest_joingroup(M))
            M.parest.hlgroupmark = line(1,1,'parent',M.grob.haflags, ...
                'marker','s','markersize',8,'linestyle','none', ...
                'color','k');
            if ~isempty(M.e)
                set(M.parest.hlgroupmark,'xdata',getflagnumber(M.e))
            else
                set(M.parest.hlgroupmark,'visible','off')
            end
            set(M.hf,'tag',''), set(M.grob.haflags,'tag','haflags') % prevent action from 'fn_imvalue'
            set(M.grob.haflags,'xtick',[],'ytick',[],'box','on')
            axis(M.grob.haflags,'image')
        end
        % actions (single vessel + group + all)
        function parest_manage(M,flag,value)
            if isempty(M.kedge), return, end
            switch flag
                % ACTIONS FOR CURRENT VESSEL
                case 'select'
                    if value>length(M.e.flow)
                        if value~=length(M.e.flow)+1 || ~M.parest.unsaved, error programming, end
                        return % no change: '(not saved)' has been selected again
                    end
                    M.parest.unsaved = false;
                    M.kflow = value;
                    M.f = M.e.flow(value); 
                    M.misc.haschanged.algo = true;
                    parest_displaypar(M)
                case 'method'
                    M.parest.unsaved = true;
                    method = M.parest.methodlist{get(M.grob.parmethod,'value')};
                    estpar = defaultpar(M,method);
                    M.f = fast_flow(M.e,method,M.f.filterpar,estpar);
                    M.misc.haschanged.algo = true;
                    parest_displaypar(M)
                case 'changefilter'
                    M.parest.unsaved = true;
                    filterpar = get(M.grob.parfilter,'string');
                    M.options.filterpar = filterpar; saveoptions(M)
                    M.f = fast_flow(M.e,M.f.method,filterpar,M.f.estpar);
                    M.misc.haschanged.algo = true;
                    parest_displaypar(M)
                case 'changepar'
                    M.parest.unsaved = true;
                    estpar = fn_str2struct(get(M.grob.parset,'string'));
                    defaultpar(M,M.f.method,estpar) % set as default
                    M.f = fast_flow(M.e,M.f.method,M.f.filterpar,estpar);
                    M.misc.haschanged.algo = true;
                    parest_displaypar(M)
                case 'active'
                    M.f.active = logical(get(M.grob.paractive,'value'));
                    M.misc.haschanged.algos = true;
                    res_display(M)
                case 'top'
                    if M.parest.unsaved, return, end
                    algolist = {M.e.flow.algo};
                    perm = [M.kflow setdiff(1:length(algolist),M.kflow)];
                    M.e.flow = M.e.flow(perm);
                    M.kflow = 1;
                    M.misc.haschanged.algos = true;
                    parest_displaypar(M)
                case 'add'
                    if ~M.parest.unsaved, return, end
                    M.parest.unsaved = false; 
                    % flow already present?
                    iflow = find(strcmp(M.f.algo,{M.e.flow.algo}),1);
                    if isempty(iflow), iflow = M.e.nest+1; end
                    % add the new flow
                    M.kflow = iflow;
                    M.e.flow(iflow) = M.f;
                    M.misc.haschanged.algos = true;
                    parest_displaypar(M)
                case 'remove'
                    if M.parest.unsaved, return, end
                    if ~isempty(M.f.files) && ~fn_reallydlg( ...
                            'Estimation results might exist,', ...
                            'are you sure you want to loose track of them?')
                        return
                    end
                    M.e.flow(M.kflow) = [];
                    M.misc.haschanged.algos = true;
                    parest_display(M)
                case 'replace'
                    % save new (unsaved) parameters and remove the last
                    % saved one, marked by M.kflow
                    if ~M.parest.unsaved, return, end
                    M.parest.unsaved = false;
                    if isempty(M.kflow), parest_manage(M,'add'), return, end
                    if ~isempty(M.e.flow(M.kflow).files) && ~fn_reallydlg( ...
                            'Estimation results might exist,', ...
                            'are you sure you want to loose track of them?')
                        return
                    end
                    % flow already present?
                    iflow = find(strcmp(M.f.algo,{M.e.flow.algo}),1);
                    if isempty(iflow)
                        iflow = M.kflow; % put new one at previous one location
                    else
                        M.e.flow(M.kflow) = []; % remove the previous one location
                        if iflow>M.kflow, iflow=iflow-1; end
                    end
                    % add the new flow
                    M.kflow = iflow;
                    M.e.flow(iflow) = M.f;
                    M.misc.haschanged.algos = true;
                    parest_displaypar(M)
                % ACTIONS FOR GROUP
                case {'group_add' 'group_replace'}
                    edgeflag = fn_switch(flag,'group_add','addedge','setedge');
                    parest_groupdefault(M,'setdefault')
                    edglist = find([M.S.edges.flagnumber]==M.e.flagnumber);
                    edglist = setdiff(edglist,M.kedge);
                    for i=edglist
                        parest_groupdefault(M,edgeflag,M.S.edges(i))
                    end
                % ACTIONS FOR ALL VESSELS
                case {'all_reorder_trg' 'all_reorder_rtg'}
                    for i=1:length(M.S.edges)
                        ei = M.S.edges(i);
                        methods = {ei.flow.method};
                        idx.t = find(strcmp(methods,'track'));
                        idx.r = find(strcmp(methods,'radon'));
                        idx.g = find(strcmp(methods,'gabor'));
                        seq = flag(end-2:end);
                        ord = [idx.(seq(1)) idx.(seq(2)) idx.(seq(3))];
                        if length(ord)~=length(ei.flow), error merde, end
                        ei.flow = ei.flow(ord);
                    end
                    M.misc.haschanged.algos = true;
                    parest_display(M)
                case 'all_addcurrent'
                    method    = M.f.method;
                    filterpar = M.f.filterpar;
                    estpar    = M.f.estpar;
                    for i=1:length(M.S.edges)
                        ei = M.S.edges(i);
                        addflow(ei,method,filterpar,estpar)
                    end
                    if M.parest.unsaved
                        M.parest.unsaved = false;
                        M.kflow = addflow(M.e,method,filterpar,estpar);
                        M.f = M.e.flow(M.kflow);
                        M.misc.haschanged.algo = true;
                        parest_display(M)
                    end
                otherwise
                    error('programming: unknown flag ''%s''',flag)
            end
            
        end
        % group
        function parest_editgroups(M)
            flaglist0 = M.S.flaglist;
            flaglist1 = fast_flagcolors(flaglist0);
            if isequal(flaglist1,flaglist0), return, end
            M.S.flaglist = flaglist1;
            % show colors table
            parest_displaygroups(M)
            % re-color vessels
            if M.im.doim
                for i=1:length(M.S.edges), vsl_color(M,i), end
            else
                scan_manage(M,'color')
            end
        end
        function parest_groupdefault(M,flag,edge)
            % function parest_defaultedge(M,'setedge|addedge',edge)
            % function parest_defaultedge(M,'setdefault')
            %---
            % manages the default estimation parameters to be assigned to a
            % newly created edge
            switch flag
                case {'addedge' 'setedge'}
                    if strcmp(flag,'setedge')
                        edge.flow(:) = [];
                    end
                    s = getflag(edge);
                    s = s.defaultpars;
                    for i=1:length(s)
                        addflow(edge,s(i).method,s(i).filterpar,s(i).estpar);
                    end
                case 'setdefault'
                    flow = M.e.flow;
                    k = M.e.flagnumber(M.S.flaglist.k);
                    s = struct( ...
                        'method',       {flow.method}, ...
                        'filterpar',    {flow.filterpar}, ...
                        'estpar',       {flow.estpar});
                    list = getflaglist(M.S);
                    list(k).defaultpars = s;
                    setflaglist(M.S,list)
                    saveoptions(M)
            end
        end
        function parest_joingroup(M,iedges,k)
            if nargin<2
                if M.im.doim
                    iedges = M.kedge; 
                else
                    iedges = M.scan.kedges; 
                end
            end
            if nargin<3
                p = get(M.grob.haflags,'currentpoint');
                k = round(p(1));
            end
            if ~isscalar(iedges), error notimplementedyet, end
            setflagnumber(M.S.edges(iedges),k)
            if M.im.doim
                vsl_color(M,iedges)
            else
                scan_manage(M,'color',iedges)
            end
            for i=iedges
                parest_groupdefault(M,'addedge',M.S.edges(i))
            end
            if ismember(M.kedge,iedges)
                set(M.parest.hlgroupmark,'xdata',getflagnumber(M.e))
                parest_display(M)
            end
        end
        % all vessels
        function parest_changeallvessels(M)
            % Prompt dialog
            % (which method)
            Prompt = {'method','method'};
            Format(1,1).format = 'integer';
            Format(1,1).type = 'list';
            methodlist = {'track','radon','gabor'};
            Format(1,1).items = methodlist;
            Format(1,1).style = 'popupmenu';
            % (remove old parameters?)
            Prompt(2,:) = {'remove old parameters','doremove'};
            Format(2,1).type = 'check';
            % (parameters to overwrite)
            Prompt(3,:) = {'parameters to overwrite','estparupd'};
            Format(3,1).type = 'edit';
            Format(3,1).size = [-1 -1];
            Format(3,1).limits = [0 10];
            % (defaults)
            DefAns = struct('method',1,'doremove',false, ...
                'estparupd','x.fieldname = value;');
            % (go)
            answer = inputsdlg(Prompt,'Parameter settings',Format,DefAns);
            if isequal(answer,DefAns), return, end % Cancel was hit
            method    = methodlist{answer.method};
            doremove  = answer.doremove;
            estparupd = fn_str2struct(answer.estparupd);
            doremoveresults = [];
            
            % Update for all vessels
            for i=1:length(M.S.edges)
                ei = M.S.edges(i);
                % (which flow exist for this method
                flow = ei.flow;
                nflow = length(flow);
                ff = strcmp({flow.method},method); % which indices to change
                % (update parameter sets)
                for k=find(ff)
                    estpark = flow(k).estpar;
                    estpark = fn_structmerge(estpark,estparupd,'skip');
                    knew = addflow(ei,method,flow(k).filterpar,estpark);
                    if knew<=nflow
                        % algorithm already existed: no new algo was
                        % created, but do not remove the old one
                        ff(knew) = false;
                    end
                end
                % (remove old parameter sets)
                if doremove
                    % (check files)
                    for k=find(ff)
                        d = dir([flow(k).savedir flow(k).filebase '_' flow(k).algo '*']);
                        if ~isempty(d) 
                            if isempty(doremoveresults)
                                answer = questdlg( ...
                                    {'Estimation results might exist,' ...
                                    'are you sure you want to loose track of them?'}, ...
                                    'warning','Yes to all','No to all','No to all');
                                doremoveresults = strcmp(answer,'Yes to all');
                            end
                            if ~doremoveresults
                                % cancel removing
                                ff(k) = false;
                            end
                        end
                    end
                    % (that's it, remove)
                    ei.flow(ff) = [];
                end
            end
            
            % Update display
            parest_displaypar(M)
        end
    end
    
    % Sections
    methods
        function sections_init(M)
            M.sections = struct;
            M.sections.active = false;
            M.sections.unsaved = [];
            M.sections.controls = [M.grob.info(5) M.grob.sectionlist ...
                M.grob.sectionlisttop M.grob.sectionlistrm M.grob.sectionlistsave M.grob.sectionlistreplace ...
                M.grob.sectionid M.grob.sectionpar M.grob.sectionestpar ...
                M.grob.sectionactive M.grob.sectionestimate M.grob.sectionestimateall];
            set(M.sections.controls,'visible','off')
        end
        function sections_menu(M)
            m = uimenu(M.hf,'label','Sections');
            item.toggle = uimenu(m,'label','show sections', ...
                'callback',@togglesections);
            uimenu(m,'label','add current to all vessels','separator','on', ...
                'callback',@(u,evt)sections_manage(M,'all_addcurrent'))
            uimenu(m,'label','replace by current in all vessels', ...
                'callback',@(u,evt)sections_manage(M,'all_replacebycurrent'))
            M.menus.sections = m;
            function togglesections(u,evt)
                b = ~M.sections.active;
                M.sections.active = b;
                set(item.toggle,'checked',fn_switch(b))
                set(M.sections.controls,'visible',fn_switch(b))
                set([M.parest.controls get(M.grob.haflags,'children')'], ...
                    'visible',fn_switch(~b))
            end
        end
        % display
        function sections_cleanup(M)
            % no selected algorithm
            M.ksection = [];
            M.u = fast_section.empty;
            % empty algo list
            set(M.grob.sectionlist,'string','','value',0)
        end
        function sections_display(M)
            if isempty(M.kedge), error programming, end
            nu = length(M.e.section);
            % any existing flow? if not, create one
            if nu>0
                M.sections.unsaved = false;
                if isempty(M.ksection)
                    M.ksection = 1;
                else
                    M.ksection = min(M.ksection,nu);
                end
                M.u = M.e.section(M.ksection); 
            else
                M.sections.unsaved = true;
                sectionpar = defaultpar(M,'section');
                M.u = fast_section(M.e,sectionpar);
            end
            drawnow
            % then, display everything
            sections_displaypar(M)
        end
        function sections_displaypar(M)
            % list of algos present in current edge
            if isempty(M.e.section)
                % MATLAB BUG!!!
                algolist = {};
            else
                algolist = {M.e.section.algo};
            end
            %disp(algolist)
            if M.sections.unsaved
                algolist{end+1} = '(not saved)';
                klist = length(algolist);
            else
                klist = M.ksection;
            end
            set(M.grob.sectionlist,'string',algolist,'value',klist)
            % ID
            set(M.grob.sectionid,'string',M.u.id)
            % section parameters
            set(M.grob.sectionpar,'string',fn_struct2str(M.u.sectionpar))
            %             % estimation parameters
            %             set(M.grob.parset,'string',fn_struct2str(M.f.estpar))
            % active flag
            set(M.grob.sectionactive,'value',M.u.active)
            % update result display
            res_display(M)
       end
       % action
       function sections_manage(M,flag,value)
            if isempty(M.kedge), return, end
            switch flag
                % ACTIONS FOR CURRENT VESSEL
                case 'select'
                    if value>length(M.e.section)
                        if value~=length(M.e.section)+1 || ~M.sections.unsaved, error programming, end
                        return % no change: '(not saved)' has been selected again
                    end
                    M.sections.unsaved = false;
                    M.ksection = value;
                    M.u = M.e.section(value); 
                    M.misc.haschanged.section = true;
                    sections_displaypar(M)
                case 'changesectionpar'
                    M.sections.unsaved = true;
                    sectionpar = fn_str2struct(get(M.grob.sectionpar,'string'));
                    defaultpar(M,'section',sectionpar) % set as default
                    try
                        M.u = fast_section(M.e,sectionpar);
                    catch
                        errordlg('failed to create section, probably vessel is too close to vessel side')
                        return
                    end
                    M.misc.haschanged.section = true;
                    sections_displaypar(M)
                case 'changeestpar'
                    %                     M.sections.unsaved = true;
                    %                     estpar = fn_str2struct(get(M.grob.sectionestpar,'string'));
                    %                     defaultpar(M,M.u.method,estpar) % set as default
                    %                     M.u = fast_flow(M.e,M.u.method,M.u.filterpar,estpar);
                    %                     M.misc.haschanged.section = true;
                    %                     sections_displaypar(M)
                case 'active'
                    M.u.active = logical(get(M.grob.paractive,'value'));
                    M.misc.haschanged.sections = true;
                    res_display(M)
                case 'top'
                    if M.sections.unsaved, return, end
                    algolist = {M.e.section.algo};
                    perm = [M.ksection setdiff(1:length(algolist),M.ksection)];
                    M.e.section = M.e.section(perm);
                    M.ksection = 1;
                    M.misc.haschanged.sections = true;
                    sections_displaypar(M)
                case 'add'
                    if ~M.sections.unsaved, return, end
                    M.sections.unsaved = false; 
                    % section already present?
                    isection = find(strcmp(M.u.algo,{M.e.section.algo}),1);
                    if isempty(isection), isection = length(M.e.section)+1; end
                    % add the new section
                    M.ksection = isection;
                    M.e.section(isection) = M.u;
                    M.misc.haschanged.sections = true;
                    sections_displaypar(M)
                case 'remove'
                    if M.sections.unsaved, return, end
                    if ~isempty(M.u.files) && ~fn_reallydlg( ...
                            'Estimation results might exist,', ...
                            'are you sure you want to loose track of them?')
                        return
                    end
                    M.e.section(M.ksection) = [];
                    M.misc.haschanged.sections = true;
                    sections_display(M)
                case 'replace'
                    % save new (unsaved) parameters and remove the last
                    % saved one, marked by M.ksection
                    if ~M.sections.unsaved, return, end
                    M.sections.unsaved = false;
                    if isempty(M.ksection), parest_manage(M,'add'), return, end
                    if ~isempty(M.e.section(M.ksection).files) && ~fn_reallydlg( ...
                            'Estimation results might exist,', ...
                            'are you sure you want to loose track of them?')
                        return
                    end
                    % section already present?
                    isection = find(strcmp(M.u.algo,{M.e.section.algo}),1);
                    if isempty(isection)
                        isection = M.ksection; % put new one at previous one location
                    else
                        M.e.section(M.ksection) = []; % remove the previous one location
                        if isection>M.ksection, isection=isection-1; end
                    end
                    % add the new section
                    M.ksection = isection;
                    M.e.section(isection) = M.u;
                    M.misc.haschanged.sections = true;
                    sections_displaypar(M)
                % ACTIONS FOR ALL VESSELS
                case {'all_addcurrent','all_replacebycurrent'}
                    sectionpar = M.u.sectionpar;
                    if strcmp(flag,'all_replacebycurrent')
                        allsection = [M.S.edges.section];
                        allotherpar = allsection(~strcmp({allsection.id},M.u.id));
                        if ~isempty([allotherpar.files]) && ~fn_reallydlg( ...
                                'Estimation results might exist,', ...
                                'are you sure you want to loose track of them?');
                            return
                        end
                        dodelete = true;
                    else
                        dodelete = false;
                    end
                    for i=1:length(M.S.edges)
                        ei = M.S.edges(i);
                        nold = length(ei.section);
                        try
                            isection = addsection(ei,sectionpar);
                        catch
                            errordlg(['failed to create section for vessel ' ...
                                num2str(i) ', probably to close to an edge of the image'])
                            ei.section = ei.section(1:nold);
                            continue
                        end
                        iother = setdiff(1:length(ei.section),isection);
                        if dodelete, ei.section(iother)=[]; end
                   end
                    if M.sections.unsaved
                        M.sections.unsaved = false;
                        M.ksection = 1;
                        M.u = M.e.section(M.ksection);
                        M.misc.haschanged.section = true;
                        M.misc.haschanged.sections = true;
                        sections_display(M)
                    end
                otherwise
                    error('programming: unknown flag ''%s''',flag)
            end
       end
    end
    
    % Trials
    methods
        function trials_init(M) %#ok<MANU>
        end
        function trials_cleanup(M)
            cla(M.grob.hatrials)
        end
        function trials_display(M)
            ha = M.grob.hatrials;
            % estimation parameters
            ISO = M.S.ISO;
            s = size(ISO);
            if s(1)==M.S.nt
                % one value per frame
                ISO = reshape(permute(ISO,[1 3 2]),[M.S.nt*M.S.nexp s(2)]);
                tidx = (1:M.S.nt*M.S.nexp)/M.S.nt + .5;
            else
                % one value per trial
                ISO = reshape(permute(ISO,[3 1 2]),[M.S.nexp s(1)*s(2)]);
                tidx = 1:M.S.nexp;
            end
            if size(ISO,2)>3
                ISO = mean(abs(ISO),2);
            end
            ymin = min(ISO(:)); ymax = max(ISO(:)); 
            m = (ymin+ymax)/2; M2 = ymax-m;
            % trial conditions
            cols = [1 1 1; .7 .7 1; 1 .7 .7; 1 .6 .6; 1 .5 .5];
            b = M.S.cond+2; % code: 1=not selected, 2=rest, 3=stim1, ...
            b(~M.S.trials)=1;
            x = shiftdim(cols(b,:),-1);
            % display
            M.trials.him = image([1 M.S.nexp],m+[-M2 M2]/2,x,'parent',ha); 
            set(ha,'ydir','normal')
            axis(ha,'normal') % cancel the 'axis image' from function fn_imvalue
            line(tidx,ISO,'parent',ha,'hittest','on', ...
                'buttondownfcn',@(h,evnt)trials_toggle(M));
        end
        function trials_setcurrent(M,itrial)
            if itrial==0, return, end
            M.ktrial = itrial;
            set(M.misc.hutrial,'value',itrial)
            if M.im.dotrialframe && isfield(M.im,'him') && any(ishandle(M.im.him))
                im_displayframe(M)
            end
            M.misc.haschanged.ktrial = true;
            res_display(M)
        end
        function trials_toggle(M)
            ha = gca;
            p = get(ha,'currentpoint');
            itrial = round(p(1));
            if strcmp(get(M.hf,'selectiontype'),'alt')
                bval = ~M.S.trials(itrial);
                b0 = M.S.cond+2; % code: 1=not selected, 2=rest, 3=stim1, ...
                fn_buttonmotion(@trial_toggle_sub)
                M.misc.haschanged.trials = true;
                res_display(M)
            else
                trials_setcurrent(M,itrial);
            end
            function trial_toggle_sub
                p = get(ha,'currentpoint');
                itrial2 = min(M.S.nexp,max(1,round(p(1))));
                idx = sort([itrial itrial2]);
                M.S.trials(idx(1):idx(2)) = bval;
                b = b0;
                b(~M.S.trials) = 1;
                cols = [1 1 1; .7 .7 1; 1 .7 .7; 1 .6 .6; 1 .5 .5];
                x = shiftdim(cols(b,:),-1);
                set(M.trials.him,'cdata',x)
            end
        end
    end
    
    % Results
    methods
        % first initialization
        function res_init(M)
            % structure
            M.res = struct( ... % structure with as many elements as graphs
                'ha',           num2cell(M.grob.hares), ...
                's',            num2cell(M.options.resdisp), ... % original structure describing what to display
                'resdisp',      [], ...     % processed structure
                'dependency',   [], ...
                'isimage',      false);     % information for display update
            % note that M.resd is initialized in res_menu because it is needed in res_menu
            % (compute dependencies)
            for k=1:7
                [dum1 dum2 dum3 M.res(k).dependency] = ...
                    res_processstructure(M,M.res(k).s);
            end
            % axes and buttons
            btpos = [-15 -15 15 15];
            for k=1:7
                axpos = get(M.res(k).ha,'position');
                uicontrol('parent',M.hf, ...
                    'userdata',k, ...
                    'string','D', ...
                    'callback',@(u,evt)res_userstructure(M,get(u,'userdata')), ...
                    'position',[axpos(1)+axpos(3) axpos(2)+axpos(4) 0 0]+btpos);
            end
        end
        function res_menu(M)
            if isempty(M.resd)
                M.resd = struct('volumescale',1,'showlegend',true,'userdisplay',false);
            end
            m = uimenu(M.hf,'label','Results');
            uimenu(m,'label','precompute results', ...
                'callback',@(u,evt)precomp(M.S,M.ktrial))
            uimenu(m,'label','clear precomputed results', ...
                'callback',@(u,evt)clearprecomp())
            uimenu(m,'label','volume scale...','separator','on', ...
                'callback',@(u,evt)setvolscal());
            item.legend = uimenu(m,'label','show legend','checked',fn_switch(M.resd.showlegend),'separator','on', ...
                'callback',@(u,evt)toggleshowlegend());
            item.user = uimenu(m,'label','enable user display','checked',fn_switch(M.resd.userdisplay),'separator','on', ...
                'callback',@(u,evt)toggleuserdisplay());
            uimenu(m,'label','edit display code...', ...
                'callback','edit fast_customresult')
            uimenu(m,'label','edit result calc. code...','separator','on', ...
                'callback',@(u,evt)opentoline(which('fast_data'),252));
            uimenu(m,'label','go back to default display','separator','on', ...
                'callback',@(u,evt)res_defaultpar(M,'getdefault'))
            uimenu(m,'label','make current display default', ...
                'callback',@(u,evt)res_defaultpar(M,'setdefault'))
            M.menus.res = m;
            function clearprecomp
                hidetmpdata(M.S.edges,'remove')
                res_display(M,[],true)
            end
            function setvolscal()
                x = fn_input('volumescaling',M.resd.volumescale,1e-2,1e4);
                if isempty(x), return, end
                M.resd.volumescale = x;
                res_display(M,[],true)
            end
            function toggleshowlegend()
                b = ~M.resd.showlegend;
                M.resd.showlegend= b;
                set(item.legend,'checked',fn_switch(b))
                res_display(M,[],true)
            end
            function toggleuserdisplay()
                b = ~M.resd.userdisplay;
                M.resd.userdisplay = b;
                set(item.user,'checked',fn_switch(b))
                res_display(M,[],true)
            end
        end
        % display
        function res_defaultpar(M,flag)
            switch flag
                case 'getdefault'
                    defresdisp = num2cell(M.options.resdisp);
                    [M.res.s] = deal(defresdisp{:});
                    res_display(M,[],true)
                case 'setdefault'
                    M.options.resdisp = [M.res.s];
                    saveoptions(M)
            end
        end
        function res_userstructure(M,i)
            M.res(i).s = fast_resdispselect(M.S,M.res(i).s);
            res_display(M,i,true)
        end
        function [x isimage isunique dependency] = res_processstructure(M,s)
            % process the 'resdisp' structure given by fast_resdispselect
            % and returns a cell array with the requested results
            % 
            % the function first produces another structure, possibly with
            % multiple elements, which will make it more easy to get the
            % result
            
            % make all field values cell arrays
            F = fieldnames(s);
            nF = length(F);
            for i=1:nF
                if ~iscell(s.(F{i})), s.(F{i})={s.(F{i})}; end
            end
            
            % initialize the intermediary structure
            x = struct( ...
                'kedge',[],'edge',[],'edgelegend',[], ...
                'data',[],'dataobject',[],'datalegend',[], ...
                'trial',[],'triallegend',[], ...
                'points',[],'pointslegend',[], ...
                'fzsub',[],'fzsublegend',[], ...
                'smooth',[],'smoothlegend',[]);
            isunique = struct('edge',true,'data',true, ...
                'trial',true,'res',true,'oddeven',true, ...
                'points',true,'fzsub',true,'smooth',true);
            dependency = struct('ktrial',false,'trials',false, ...
                'kedge',false,'algo',false,'algos',false,'fastflow',false, ...
                'section',false,'sections',false);
            isimage = false; 
            
            % Edge: put the adequate fast_edge object
            kedges = s.edge;
            ne = length(s.edge);
            for i=1:ne
                if isequal(kedges{i},'current')
                    dependency.kedge = true;
                    if isempty(M.kedge), x=[]; return, end % no edge
                    kedges{i}=M.kedge; 
                else
                    % check that edge exist
                    if ~isnumeric(kedges{i}), error('wrong edge'), end
                    if kedges{i}(1)<0
                        % groupe numbers!
                        kgroup = -kedges{i}; kedges{i} = [];
                        flagnumbers = getflagnumber(M.S.edges);
                        flagnumbers(~[M.S.edges.active])=0;
                        for j=1:length(kgroup)
                            kedges{i} = [kedges{i} find(flagnumbers==kgroup(j))];
                        end
                    end
                    if ~all(ismember(kedges{i},1:length(M.S.edges)))
                        error('wrong edge')
                    end
                end
            end
            kedges = unique([kedges{:}]); % note that this works also when M.kedge is empty!
            ne = length(kedges);
            isunique.edge = (ne==1);
            x = repmat(x,1,ne);
            for i=1:ne
                x(i).kedge = kedges(i);
                x(i).edge = M.S.edges(kedges(i));
                x(i).edgelegend = ['vessel ' num2str(kedges(i))];
            end
            
            % Data/Algo: that is difficult!!!
            ndata = length(s.data);
            if ndata>1, isunique.data = false; end
            xdata = cell(ndata,ne); % the number of cases for each data flag will depend on this flag and on which vessel
            % loop on edges 
            for j=1:ne
                % loop on data (set data)
                for i=1:ndata
                    y = {};
                    z0 = x(j); 
                    if strcmp(s.data{i},'width'), z0.data='sigma'; else z0.data=s.data{i}; end
                    if strcmp(s.data{i},'result'), nmethod=length(s.algo); else nmethod=1; end
                    if nmethod>1, isunique.data = false; end
                   % loop on method (set dataobject and datalegend)
                    for k=1:nmethod
                        z = z0;
                        method = s.algo{k};
                        if strcmp(s.data{i},'volume') || strcmp(s.data{i},'flow')
                            z.dataobject = z.edge;
                            z.datalegend = s.data{i};
                            y{end+1} = z; %#ok<AGROW>
                        elseif strcmp(s.data{i},'section') || strcmp(s.data{i},'width')
                            dependency.section = true;
                            if x(j).kedge==M.kedge
                                z.dataobject = M.u;
                            else
                                if isempty(x(j).edge.section), break, end % no result to display for this vessel
                                z.dataobject = x(j).edge.section(1);
                            end
                            z.datalegend = z.dataobject.algo;
                            if strcmp(s.data{i},'width'), z.datalegend = [z.datalegend '-width']; end
                            y{end+1} = z; %#ok<AGROW>
                        elseif strcmp(method,'current')
                            dependency.algo = true;
                            if x(j).kedge==M.kedge
                                z.dataobject = M.f;
                            else
                                if isempty(x(j).edge.flow), break, end % no result to display for this vessel
                                z.dataobject = x(j).edge.flow(1);
                            end
                            z.datalegend = z.dataobject.algo;
                            y{end+1} = z; %#ok<AGROW>
                        else
                            dependency.algos = true;
                            doall = any(findstr(method,'all'));
                            if doall, method = strrep(method,'all',''); end
                            if isempty(x(j).edge.flow), continue, end % no result to display for this vessel and this method
                            iflow = find(strcmp(method,{x(j).edge.flow.method}));
                            if isempty(iflow), continue, end % no result to display for this vessel and this method
                            if ~doall, iflow = iflow(1); end
                            if length(iflow)>1, isunique.data = false; end
                            for l=1:length(iflow)
                                fl = x(j).edge.flow(iflow(l));
                                if ~fl.active, continue, end
                                z.dataobject = fl;
                                if strcmp(s.data{i},'volumef'), z.datalegend='volumef'; else z.datalegend = fl.algo; end
                                y{end+1} = z; %#ok<AGROW>
                            end
                        end
                    end
                    if isempty(y), continue, end
                    y = [y{:}];
                    xdata{i,j} = y;
                end
            end
            x = [xdata{:}];
            if isempty(x), return, end
            if isunique.data && ~isscalar(x) && strcmp(x(1).data,'result')
                % several result displays: same algo?
                datalegends = strvcat(x.datalegend); %#ok<VCAT>
                difflegends = [any(diff(datalegends,1,1),1) 1];
                idx = find(difflegends,1,'first'); % first character where legends differ
                datalegend = x(1).datalegend(1:idx-1);
                if isempty(datalegend), datalegend = 'result'; end
                [x.datalegend] = deal(datalegend);
            end
                        
            % Trials: replace 'current' by M.ktrial, duplicate x if
            % necessary
            y = x; x = {};
            if length(s.trial)>1, isunique.trial = false; end
            idx = 0;
            for i=1:length(s.trial)
                trialflag = s.trial{i};
                if isequal(trialflag,'current')
                    trialflag=M.ktrial;
                    dependency.ktrial = true;
                end
                % single trial or average over time -> easy
                if isnumeric(trialflag) || strcmp(trialflag,'timeavg')
                    [y.trial] = deal(trialflag);
                    if isnumeric(trialflag)
                        [y.triallegend] = deal(['trial ' num2str(trialflag)]);
                    else
                        [y.triallegend] = deal('all experiment');
                    end
                    idx = idx+1; 
                    x{idx} = y; %#ok<AGROW>
                    continue
                end
                % average -> need to find which trials to average
                if ~strcmp(trialflag,'average'), error('wrong trial flag ''%s''',trialflag), end
                [y.trial] = deal('average');
                nres = length(s.res); noddeven = length(s.oddeven); 
                if nres>1,     isunique.res = false;     end
                if noddeven>1, isunique.oddeven = false; end
                isunique.trial = isunique.res & isunique.oddeven;
                for k=1:nres
                    for j=1:noddeven
                        result = s.res{k}; oddeven = s.oddeven{j}; 
                        [y.res]     = deal(result);
                        [y.oddeven] = deal(oddeven);
                        if strcmp(result,'all')
                            if strcmp(oddeven,'all')
                                [y.triallegend] = deal('all trials');
                            else
                                [y.triallegend] = deal(oddeven);
                            end
                        else
                            if strcmp(oddeven,'all')
                                [y.triallegend] = deal(result);
                            else
                                [y.triallegend] = deal([result '-' oddeven]);
                            end
                        end
                        idx = idx+1;
                        x{idx} = y; %#ok<AGROW>
                    end
                end
                dependency.trials = true;
            end
            x = [x{:}];
            
            % Points
            isimage = deal(isscalar(s.points) && strcmp(s.points{1},'image'));
            npoints = length(s.points);
            isunique.points = (npoints==1);
            x = repmat(x,npoints,1);
            for i=1:npoints
                [x(i,:).points]       = deal(s.points{i});
                [x(i,:).pointslegend] = deal(s.points{i});
            end
            x = x(:)';
            
            % Smoothing, Frame zero subtraction
            G = {'smooth' 'fzsub'};
            for k=1:length(G)
                g = G{k};
                ng = length(s.(g));
                isunique.(g) = (ng==1);
                x = repmat(x,ng,1);
                for i=1:ng
                    [x(i,:).(g)]            = deal(s.(g){i});
                    [x(i,:).([g 'legend'])] = deal(num2str(s.(g){i}));
                end
                x = x(:)';
            end
            
            % Special - Volume scaling
            if M.resd.volumescale
                [x.volumescale] = deal(M.resd.volumescale);
            end
            
            % Special - Width for image display of vessel section
            if isimage && ~isscalar(x), error programming, end
            if isimage && strcmp(x.data,'section')
                x(2) = x;
                x(2).data = 'profile';
            end
            
            % Special - Color for display
            nx = length(x);
            if isunique.oddeven
                for i=1:nx, x(i).kcolor=i; end
            else
                oddeven = {x.oddeven};
                kodd  = find(strcmp(oddeven,'odd'));
                keven = find(strcmp(oddeven,'even'));
                for i=1:nx/2, x(kodd(i)).kcolor=i; x(keven(i)).kcolor=i; end
            end
            
            % ouf!!!
        end
        function res_display(M,igraph,forceflag)
            global Y
            % display holding?
            if M.misc.holdresdisplay, return, end
            % input
            doall = nargin<2 || isempty(igraph);
            if doall, igraph = 1:length(M.res); end
            if nargin<3, forceflag = false; end
            % update the display specifications if necessary
            if M.misc.haschanged.fastflow
                for i=igraph % igraph should be 1:length(M.res)
                    M.res(i).s = fast_resdispselect(M.S,M.res(i).s,'correction');
                end
                forceflag = true;
            end
            % Loop on displays
            for i=igraph
                resi = M.res(i);
                % should display be updates?
                depchg = [resi.dependency M.misc.haschanged];
                depchg = permute(cell2mat(struct2cell(depchg)),[3 2 1]); % 2 x nfields array of logicals
                doupd = forceflag || any(all(depchg));
                if ~doupd, continue, end
                % process the graph structure first
                [resdisp M.res(i).isimage isunique M.res(i).dependency] = ...
                    res_processstructure(M,resi.s);
                % get the data
                nx = length(resdisp);
                x = cell(1,nx);
                ok = false(1,nx);
                for k=1:nx
                    objk = resdisp(k).dataobject;
                    x{k} = getresult(objk,resdisp(k)); % dataobject can be a fast_edge or fast_flow object
                    ok(k) = ~isempty(x{k});
                    % special additional computation
                    if ~ok(k) && ~isempty(findstr(resdisp(k).data,'volume')) ...
                            && isequal(resdisp(k).trial,M.ktrial)
                        % TD 07/10/2015 - previous code below seems buggy
                        %                         edgek = fn_switch(resdisp(k).data,'volume',objk, ...
                        %                             'volumef',objk.edge);
                        switch class(objk)
                            case 'fast_edge'
                                edgek = objk;
                            case 'fast_flow'
                                edgek = objk.edge;
                        end
                        % end TD 07/10/2015
                        % no volume data -> try to get it by interpolating the movie data
                        switch edgek.interpflag
                            case 'movie'
                                success = loadtrial(M.S,M.ktrial);
                                if ~success, continue, end
                                I = fast_interptube(Y,edgek);
                            case 'scan'
                                fname = [edgek.S.datadir deblank(edgek.filepar.files(M.ktrial,:))];
                                I = fn_readimg(fname);
                            case 'fake'
                                I = fast_fakedata(edgek.fakepar);
                            otherwise
                                error programming
                        end
                        setdatatmp(edgek,M.ktrial,I)
                        x{k} = getresult(objk,resdisp(k)); % dataobject can be a fast_edge or fast_flow object
                        ok(k) = ~isempty(x{k});
                        if ~ok(k), error programming, end
                    elseif ~ok(k) && strcmp(resdisp(k).data,'section') ...
                            && isequal(resdisp(k).trial,M.ktrial)
                        % no section data -> try to get it by interpolating the movie data
                        if ~M.im.doim, error programming, end
                        success = loadtrial(M.S,M.ktrial);
                        if ~success, continue, end
                        U = fast_vesselsection(Y,objk.edge,objk.sectionpar);
                        setsectiondatatmp(objk,M.ktrial,U)
                        x{k} = getresult(resdisp(k).dataobject,resdisp(k)); % dataobject can be a fast_edge or fast_flow object
                        ok(k) = ~isempty(x{k});
                        if ~ok(k), error programming, end
                    elseif ok(k) && strcmp(resdisp(k).data,'result') ...
                            && isfield(resdisp(k).dataobject.user,'resultscale')
                        % scale result using user-defined 'resultscale' value
                        x{k} = x{k}*resdisp(k).dataobject.user.resultscale;
                    end
                end
                nx = sum(ok);
                resdisp = resdisp(ok);
                M.res(i).resdisp = resdisp;
                % cell array -> array(s)
                if M.res(i).isimage
                    if ~isscalar(x)
                        if length(x)~=2, error programming, end
                        y = x{2};
                    else
                        y = {};
                    end
                    x = x{1};
                else
                    try
                        x = [x{:}];
                    catch %#ok<*CTCH>
                        % make all time courses the same size (correct for
                        % binning, all x{k} assumed to be vectors)
                        x = x(ok); nx = length(x);
                        nt = zeros(1,nx);
                        for k=1:nx, nt(k)=length(x{k}); end
                        bin = max(nt)./nt;
                        for k=1:nx, x{k}=kron(x{k},ones(bin(k),1)); end
                        x = [x{:}];
                    end % if sizes are different, x remains a cell array
                end
                % clear display?
                if isempty(x)
                    cla(resi.ha), legend(resi.ha,'off'), title(resi.ha,'')
                    continue
                end
                % display
                if M.res(i).isimage
                    % IMAGE DISPLAY
                    % shift pixels?
                    if ismember(resdisp(1).data,{'volume','volumef','lines'})
                        npix = get(M.grob.hulinshift,'value')-1;
                        if npix>10, npix=-(npix-10); end
                        if npix
                            [np nt] = size(x);
                            np2 = np+npix;
                            [ii jj] = ndgrid(1:np,1:nt);
                            ii = fn_add(ii,-(0:nt-1)*npix); % this is the pixel shift
                            ii = 1+mod(ii-1,np2);
                            idx = ii + np2*(jj-1);
                            tmp = ones(np2,nt)*min(x(:));
                            tmp(idx) = x;
                            x = tmp;
                        end
                    end
                    % display (take into account previous zooming)
                    ax = axis(resi.ha);
                    imagesc(x,'parent',resi.ha)
                    if ax(1), try fn_imvalue('chgx',ax(1:2),resi.ha), end, end %#ok<TRYNC>
                    % if strcmp(resdisp(1).data,'section'), axis(resi.ha,'normal'), else axis(resi.ha,'image'), end
                    % title
                    title(resi.ha,[resdisp(1).datalegend ' (' ...
                        resdisp(1).edgelegend ', ' resdisp(1).triallegend ')'], ...
                        'fontsize',10)
                    % width of the vessel section
                    if strcmp(resdisp(1).data,'section') && ~isempty(y)
                        halfwidth = y(1,:)*1.18; % sigma -> HWHM
                        center = y(2,:);
                        tt = 1:size(x,2);
                        hl=line([tt NaN tt],[center+halfwidth NaN center-halfwidth], ...
                            'color','b','parent',resi.ha);
                        set(hl,'hittest','on','buttondownfcn',@(hl,e)delete(hl))
                    end
                else
                    % GRAPH DISPLAY
                    hl = plot(x,'parent',resi.ha);
                    % color
                    colors = 'bgrcmykw';
                    ncol = length(colors);
                    for k=1:nx
                        col = colors(1+mod(resdisp(k).kcolor-1,ncol));
                        set(hl(k),'color',col)
                    end
                    % enable trial selection by clicking on the curve?
                    if isequal(resdisp(1).trial,'timeavg')
                        set(hl,'hittest','on','buttondownfcn',@(u,evt)trials_toggle(M))
                    end
                    % make title and legend!
                    nres = length(resdisp);
                    leg  = cell(1,nres);
                    sepleg = '';
                    if isunique.data
                        name1 = resdisp(1).datalegend;
                    else
                        name1 = 'compare';
                        for k=1:nres, leg{k} = resdisp(k).datalegend; end
                        sepleg = ', ';
                    end
                    if isunique.edge
                        name2 = resdisp(1).edgelegend;
                    else
                        nvessel = length(unique([resdisp.kedge]));
                        name2 = [num2str(nvessel) ' vessels'];
                        for k=1:nres, leg{k} = [leg{k} sepleg resdisp(k).edgelegend]; end
                        sepleg = ', '; 
                    end
                    if isunique.trial
                        name3 = resdisp(1).triallegend;
                    else
                        if isunique.res
                            name3 = resdisp(1).res;
                            if strcmp(name3,'all'), name3 = 'all trials'; end
                        else
                            name3 = 'rest+stim';
                        end
                        for k=1:nres, leg{k} = [leg{k} sepleg resdisp(k).triallegend]; end
                        sepleg = ', ';
                    end
                    if ~isunique.points
                        for k=1:nres, leg{k} = [leg{k} sepleg resdisp(k).pointslegend]; end
                        sepleg = ', ';
                    end
                    if isunique.smooth
                        if resdisp(1).smooth, name3 = [name3 ', smooth']; end %#ok<AGROW>
                    else
                        for k=1:nres, leg{k} = [leg{k} sepleg resdisp(k).smoothlegend]; end
                    end
                    if isunique.fzsub
                        if resdisp(1).fzsub, name3 = [name3 ', fzsub']; end %#ok<AGROW>
                    else
                        for k=1:nres, leg{k} = [leg{k} sepleg resdisp(k).fzsublegend]; end
                        sepleg = ', ';
                    end
                    title(resi.ha,[name1 ' (' name2 ', ' name3 ')'],'fontsize',10)
                    if M.resd.showlegend && ~isempty(sepleg)
                        legend(hl,leg,'Location','NorthWest')
                    end
                end
            end
            F = fieldnames(M.misc.haschanged);
            for i=1:length(F), M.misc.haschanged.(F{i})=false; end
            % User-defined display
            if M.resd.userdisplay && doall
                fast_customresult
            end
        end
        % actions
        function res_trajectories(M,ha)
            p = get(ha,'currentpoint'); p = p(1,1:2)';
            % list of estimation results
            flowlist = M.e.flow;
            if M.parest.unsaved, flowlist(end+1) = M.f; end
            nx = M.f.np; nt = M.f.nt;
            flowlist = flowlist([flowlist.active] ...
                & [flowlist.np]==nx & [flowlist.nt]==nt);
            imlist = find([M.res.isimage]);
            % list of colors
            colors = 'bgrcmykw';
            ncol = length(colors);
            % size parameters
            npix = get(M.grob.hulinshift,'value')-1;
            if npix>10, npix=-(npix-10); end
            np2 = nx+npix;
            % compute trajectories
            for i_res=1:length(flowlist)
                V = getresult(flowlist(i_res),M.ktrial);
                if isempty(V), continue, end
                t = round(p(1)); x = p(2);
                if npix 
                    % pix shift!
                    x = 1+mod((x-1)+t*npix,np2); x = min(x,nx); 
                end
                t0 = t; x0 = x; xa = []; xb = [];
                while t>0 && x>=1 && x<=nx
                    vxt = interp1(1:nx,V(:,t),x);
                    x = x - vxt;
                    t = t-1;
                    xa(end+1) = x; %#ok<AGROW>
                end
                t = t0; x = x0;
                while t<=nt && x>=1 && x<=nx
                    vxt = interp1(1:nx,V(:,t),x);
                    x = x + vxt;
                    t = t+1;
                    xb(end+1) = x; %#ok<AGROW>
                end
                tdata = t0 + (-length(xa):length(xb));
                xdata = [xa(end:-1:1) x0 xb];
                if npix
                    % pix shift!
                    xdata = 1+mod(xdata-tdata*npix,np2);
                end
                % display
                col = colors(1+mod(i_res-1,ncol));
                hl = [];
                for i_ha=imlist
                    hl(end+1) = line(tdata,xdata, ...
                        'color',col, ...
                        'parent',M.grob.hares(i_ha)); %#ok<AGROW>
                end
                set(hl,'tag','trajectory','hittest','on', ...
                    'buttondownfcn',@(h,evnt)trydelete(hl,M.hf))
            end
        end
        function res_estimate(M,flag)
            if isempty(M.e), return, end
            doall = any(findstr(flag,'all'));
            force = ~doall;
            switch flag
                case 'flow'
                    objs = M.f;
                case 'flowall'
                    objs = M.e.flow;
                    objs = objs([objs.active]);
                case 'section'
                    objs = M.u;
                case {'allsection','sectionall'} % TODO: remove first case
                    objs = M.e.section;
                    objs = objs([objs.active]);
                otherwise 
                    error flag
            end
            flag=strrep(flag,'all',''); % 'flow' or 'section'
            if isempty(objs), return, end
            for obj = objs
                disp(obj.algo), tic
                if force || isempty(getresult(obj,M.ktrial))
                    % estimate
                    switch flag
                        case 'flow'
                            flow = obj;
                            method = flow.method;
                            estpar = flow.estpar;
                            J0V0 = {};
                            % filter
                            I = getdataf(flow,M.ktrial);
                            if isempty(I), error('problem: no data'), end
                            % check fasttrack initialization
                            if strcmp(method,'track') && ischar(estpar.v0)
                                switch estpar.v0
                                    case 'last'
                                        estpar.v0 = [];
                                    case 'radon'
                                        estpar.v0 = [];
                                        kradon = find(strcmp({M.e.flow.method},'radon'));
                                        if isempty(kradon)
                                            msgbox('do radon method')
                                            return
                                        end
                                        fradon = M.e.flow(kradon(1));
                                        V0 = getresult(fradon,M.ktrial);
                                        if isempty(V0)
                                            fprintf([repmat('\b',1,length(flow.algo+1)) ...
                                                fradon.algo '\n']), tic
                                            V = fast_radon(I,fradon.estpar);
                                            setresulttmp(fradon,M.ktrial,V);
                                            V0 = V;
                                            fprintf('\b (%is)\n%s\n',floor(toc),flow.algo), tic
                                        end
                                        J0V0 = {I V0};
                                    case {'fft','est'}
                                        % keep it
                                    otherwise
                                        error merde
                                end
                            end
                            % estimate
                            switch method
                                case 'track'
                                    [J V] = fast_track(I,estpar,J0V0{:});
                                    setresulttmp(flow,M.ktrial,V,J);
                                otherwise
                                    V = feval(['fast_' method],I,estpar);
                                    setresulttmp(flow,M.ktrial,V);
                            end
                        case 'section'
                            % section data
                            U = getsectiondata(obj,M.ktrial);
                            if isempty(U), error('problem: no section data'), end
                            % estimate
                            profile = fast_vesselwidth(U);
                            setresulttmp(obj,M.ktrial,profile)
                        otherwise
                            error programming
                    end
                    fprintf('\b (%is)\n',floor(toc))
                else
                    fprintf('\b (result exists already)\n')
                end
            end
            % display result
            switch flag
                case 'flow'
                    M.misc.haschanged.algos = true;
                    M.misc.haschanged.algo = true;
                case 'section'
                    M.misc.haschanged.sections = true;
                    M.misc.haschanged.section = true;
            end
            res_display(M)
        end
    end

    % Misc
    methods
        function par = defaultpar(M,method,value)
            % function par = defaultpar(M,method)
            % function defaultpar(M,method,value)
            %---
            % get or set default parameters for specific estimation method
            if nargin==2
                % GET
                if isfield(M.options.parest,method)
                    par = M.options.parest.(method);
                else
                    % default parameters given by the function itself
                    switch method
                        case 'track'
                            par = fast_track('par');
                        case 'radon'
                            par = fast_radon('par');
                        case 'gabor'
                            par = fast_gabor('par');
                        case 'section'
                            par = fast_vesselsection('par');
                        otherwise
                            error('unknown method ''%s''',method)
                    end
                    M.options.parest.(method) = par;
                    saveoptions(M)
                end
            else
                % SET
                M.options.parest.(method) = value;
                saveoptions(M)
            end
        end
        function misc_axesclick(M,ha)
            % check that current point is inside an axes
            p = get(ha,'currentpoint'); p = p(1,1:2)';
            ax = axis(ha);
            if any(p'<ax([1 3]) | p'>ax([2 4])), return, end % cancel if point is outside of axes
            if ha==M.grob.haim
                vessels_selectpoly(M,ha)
            elseif ismember(ha,M.grob.hares) && strcmp(get(M.hf,'selectionType'),'open')
                res_trajectories(M,ha)
            end
        end
        function access(M) %#ok<MANU>
            keyboard
        end
    end
    
    
end



%------
% TOOLS
%------

% delete trajectory lines 
function trydelete(hl,hf)
switch get(hf,'selectionType')
    case 'normal'
        delete(hl(ishandle(hl)))
    case 'alt'
        delete(findobj(hf,'tag','trajectory'))
end
end
    
function u = togglegroup(labels,fun,groupprop,buttonprop)
% function u = togglegroup(labels,fun,groupprop,buttonprop)
%---
% Input:
% - labels      cell array of char array - button labels
% - fun         function handle @fun(str) - function to be executed, should
%               have one argument (the label)
% - groupprop   cell array of option/value pairs - properties for the group
% - buttonprop  cell array of option/value pairs - properties for the
%               buttons


% input
if nargin<3, groupprop={}; end
if nargin<4, buttonprop={}; end

% container
u = uibuttongroup('selectionchangefcn', ...
    @(hu,evt)feval(fun,get(get(hu,'selectedobject'),'string')), ...
    groupprop{:});

% elements
n = length(labels);
for i=1:n
    uicontrol('style','togglebutton','parent',u, ...
        'string',labels{i}, ...
        'units','normalized','position',[(i-1)/n 0 1/n 1], ...
        buttonprop{:})
end

end

% TO TRASH
function k = getflagnumber(E)
klist = E(1).S.flaglist.k; % index of current flag list
k = cat(1,E.flagnumber);
k = k(:,klist)';
end
