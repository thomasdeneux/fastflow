classdef fast_edgesmouse < interface
    properties
        S
        parameters
        menuitems
        kY
        anchors
        vsls
        vslmenu
        selpolyflag
        tmp
        curvsl
        ha
        hb
        localparinput
        noteinput
        statusbar
        resdisp
        gabor = struct('par',[],'np',0,'FGabor',[]);
    end
    
    % Constructor
    methods
        function E = fast_edgesmouse(S)
            % init interface
            E = E@interface(S.hf+1,'VESSEL SELECTION');
            fn_imvalue image
            
            % init parameters from S
            E.parameters = S.parameters;
            
            % data from fastflow object
            E.S = S;
            
            % check data
            global Y
            if ~iscell(Y) && ~isequal(size(Y),[S.nx S.ny S.nt])
                error programming
            end
            
            % edges
            edges = S.edges;
            if all(all(any(S.ISO)))
                m = min(min(S.ISO,[],3));
                M = max(max(S.ISO,[],3));
                % qui va ou: au pif au metre!
                contour = [3-m(3) S.nx-M(3)-2 3-m(2) S.ny-M(2)-2];
            else
                contour = [];
            end
            % note that if 'fake' data, a menu will be created by
            % init_menus, called by interface_end
            E.kY = 1;
            
            % images
            E.ha(1) = axes;
            imagesc(S.CS'), colormap gray
            addlistener(E.ha(1),'YLim','PostSet',@(m,evnt)autoclip(E.ha(1),S.CS));
            
            E.ha(2) = axes;
            maskcol = E.S.mask' / (2*std(E.S.mask(:)));
            maskcol = 1+max(0,min(255,round((maskcol+1)/2*255)));
            cm = signcheck(256);
            [ni nj] = size(maskcol);
            maskcol = reshape(cm(maskcol(:),:),[ni nj 3]);
            image(maskcol)
            set(E.ha,'xtick',[],'ytick',[])
            set(E.hf,'windowbuttondownfcn',@(hf,evnt)axesclick(E,gca))
            
            % contour
            if ~isempty(contour)
                xdata = contour([1 2 2 1 1]);
                ydata = contour([3 3 4 4 3]);
                line(xdata,ydata,'parent',E.ha(1),'color','k')
                line(xdata,ydata,'parent',E.ha(2),'color','b')
            end
            
            % lines
            E.hb(1) = axes;
            E.hb(2) = axes;
            E.hb(3) = axes;
            E.hb(4) = axes;
            E.hb(5) = axes;
            E.resdisp = {'I' 'J' 'V' 'V1' 'V3'};
            
            % estimation parameters
            set(E.hf, ...
                'defaultuicontrolhorizontalalignment','left', ...
                'defaultuicontrolfontsize',9)
            hu(1)=uicontrol('style','text','string',' VESSEL DETECTION');
            hu(2)=uicontrol('style','edit','max',2, ...
                'backgroundcolor','w', ...
                'string',E.parameters.segment, ...
                'callback',@(u,evnt)chgsegmentpar(E,get(u,'string')));
            % line filtering parameters
            hu(3)=uicontrol('style','text','string',' LINES FILTERING');
            hu(4)=uicontrol('style','edit','max',2, ...
                'backgroundcolor','w', ...
                'string',E.parameters.linefilter, ...
                'callback',@(u,evnt)chglinefilter(E,get(u,'string')));
            % line tracking parameters
            hu(5)=uicontrol('style','text','string',' LINES TRACKING (GLOBAL)');
            hu(6)=uicontrol('style','edit','max',2, ...
                'backgroundcolor','w', ...
                'string',E.parameters.estpar, ...
                'callback',@(u,evnt)chglinetrackglobal(E,get(u,'string')));
            % line tracking parameters (overwrite for individual vessel)
            hu(7)=uicontrol('style','text','string',' LINES TRACKING (THIS VESSEL)');
            hu(8)=uicontrol('style','edit','max',2, ...
                'backgroundcolor','w', ...
                'string','', ...
                'callback',@(u,evnt)chglinetracklocal(E,get(u,'string')));
            % result display parameters
            hu(11)=uicontrol('style','text','string',' RESULTS DISPLAY');
            hu(12)=uicontrol('style','edit','max',2, ...
                'backgroundcolor','w', ...
                'string',E.parameters.resdisp, ...
                'callback',@(u,evnt)chgresdisppar(E,get(u,'string')));
            % help
            hu(9)=uicontrol('style','text','max',2, ...
                'string',{'HELP: - Double-click image to set a new anchor point.', ...
                '- To initialize a vsl, click an anchor point (start), then intermediary points, then and anchor point (end).', ...
                '- Vessel middle-points have context menu for more options.'});
            % info display
            hu(10)=uicontrol('style','text');
            % movie
            if iscell(Y)
                fun = @(u,evnt)fn_movie(Y{E.kY});
            else
                fun = @(u,evnt)fn_movie(Y);
            end
            hu(13)=uicontrol('string','MOVIE','callback',fun);
            
            % add a note to vessels
            hu(14) = uicontrol('style','text','fontsize',8, ...
                'string','note:');
            hu(15) = uicontrol('style','edit','max',1, ...
                'fontsize',8,'horizontalalignment','left', ...
                'backgroundcolor','w', ...
                'callback',@(hu,evnt)setnote(get(hu,'string')));
            E.noteinput = hu(15);
            
            % local function 
            function setnote(str)
                if ~isempty(E.curvsl)
                    E.vsls(E.curvsl).edge(1).note = str;
                end
            end
           
            E.localparinput = hu(8);
            E.statusbar = hu(10);
            
            % buttons
            hx(1) = uicontrol('string','SAVE', ...
                'callback',@(u,evnt)savebutton(E));
            % dne = hu(6);
            hx(2) = uicontrol('string','FRAMES', ...
                'callback',@(u,evnt)chgframepositions(E));
            hx(3) = uicontrol('string','ESTIMATE', ...
                'callback',@(u,evnt)linestrack(E));
            hx(4) = uicontrol('string','FILTER', ...
                'callback',@(u,evnt)refilter(E));
            
            % change positions in figure
            E.grob = struct('hf',E.hf,'ha',E.ha,'hb',E.hb, ...
                'hu',hu,'hx',hx);
            interface_end(E)
            
            % edge context menus
            m = uicontextmenu('parent',E.hf);
            E.vslmenu = m;
            uimenu(m,'label','flag good', ...
                'callback',@(hp,evnt)vslflag(E,get(gco,'userdata'),'good'))
            uimenu(m,'label','flag medium', ...
                'callback',@(hp,evnt)vslflag(E,get(gco,'userdata'),'medium'))
            uimenu(m,'label','flag bad', ...
                'callback',@(hp,evnt)vslflag(E,get(gco,'userdata'),'bad'))
            uimenu(m,'label','set/remove for analysis','separator','on', ...
                'callback',@(hp,evnt)toggleactive(E,get(gco,'userdata')))
            uimenu(m,'label','bad, remove from analysis', ...
                'callback',@(hp,evnt)vslbad(E,get(gco,'userdata')))
            uimenu(m,'label','upside-down','separator','on', ...
                'callback',@(hp,evnt)upsidedown(E,get(gco,'userdata')))
            uimenu(m,'label','re-segment edge', ...
                'callback',@(hp,evnt)resetvessel(E,get(gco,'userdata')))
            uimenu(m,'label','delete edge','separator','on', ...
                'callback',@(hp,evnt)vsldelete(E,get(gco,'userdata')))
            
            % graph of anchors and edges
            % (find anchors)
            ne = size(edges,2);
            anchorpos = zeros(2,2*ne);
            for i=1:ne
                anchorpos(:,2*i-1:2*i) = edges(1,i).snake(:,[1 end]);
            end
            [anchorpos m2u u2m] = unique(anchorpos','rows');
            anchorpos = anchorpos';
            u2m = reshape(u2m,2,ne);
            % (anchors structure and display)
            na = size(anchorpos,2);
            E.anchors = [];
            for i=1:na, anchornew(E,anchorpos(:,i)), end
            % (edges structure and display)
            E.vsls = [];
            for i=1:ne, vslnew(E,edges(:,i),u2m(:,i)), end
            
            % other fields in E
            E.selpolyflag = '';
            E.tmp = [];
            E.curvsl = [];
            assignin('base','E',E)
        end
        function init_menus(E)
            delete(findobj(E.hf,'type','uimenu'))
            init_menus@interface(E)
            
            % some options
            m = E.menus.interface;
            items = struct;
            items.immest = uimenu(m,'label','immediate estimate','separator','on', ...
                'checked',fn_switch(E.options.doresponsiveest), ...
                'callback',@immestswitch);
            function immestswitch(u,evt) %#ok<INUSD>
                b = ~E.options.doresponsiveest;
                E.options.doresponsiveest = b;
                saveoptions(E)
                set(items.immest,'checked',fn_switch(b))
            end
            E.menuitems.edgesmouse = items;
            
            % menu to switch between different fake data
            if status(E.S,'fake')
                % create a menu to choose which data to use
                s = E.S.fake;
                m = uimenu(E.hf,'label','fake parameters');
                items = cell(1,4);
                for i=1:s.nst
                    mi=uimenu(m,'label',['step ' num2str(s.step(i))], ...
                        'callback',@(u,evnt)chgFakeIdx(E,1,i));
                    if i==1, set(mi,'checked','on'), end
                    items{1}(end+1)=mi;
                end
                for i=1:s.nth
                    mi=uimenu(m,'label',['thresh ' num2str(s.thresh(i))], ...
                        'callback',@(u,evnt)chgFakeIdx(E,2,i));
                    if i==1, set(mi,'checked','on','separator','on'), end
                    items{2}(end+1)=mi;
                end
                for i=1:s.nji
                    mi=uimenu(m,'label',['jitter ' num2str(s.jitter(i))], ...
                        'callback',@(u,evnt)chgFakeIdx(E,3,i));
                    if i==1, set(mi,'checked','on','separator','on'), end
                    items{3}(end+1)=mi;
                end
                for i=1:s.nsn
                    mi=uimenu(m,'label',['shotnoise ' num2str(s.shotnoise(i))], ...
                        'callback',@(u,evnt)chgFakeIdx(E,4,i));
                    if i==1, set(mi,'checked','on','separator','on'), end
                    items{4}(end+1)=mi;
                end
                % possibility to re-generate data
                uimenu(m,'label','generate new data','separator','on', ...
                    'callback',@(u,evnt)chgFakeData(E,false))
                % save handles
                E.menus.fake = m;
                E.menuitems.fake = items;
            end
            
            % menus to select what to see in the results graphs
            for k = 1:5
                m = uimenu(E.hf,'label',['result display ' num2str(k)]);
                % copy content
                uimenu(m,'label','copy content', ...
                    'callback',@(u,ebnt)copyDisplay(E.grob.hb(k)))
                % choose what to display
                s.I=uimenu(m,'label','I','separator','on', ...
                    'callback',@(u,evnt)chgResDisp(E,k,'I'));
                s.J=uimenu(m,'label','J', ...
                    'callback',@(u,evnt)chgResDisp(E,k,'J'));
                s.V=uimenu(m,'label','V','separator','on', ...
                    'callback',@(u,evnt)chgResDisp(E,k,'V'));
                s.G=uimenu(m,'label','G', ...
                    'callback',@(u,evnt)chgResDisp(E,k,'G'));
                s.R=uimenu(m,'label','R', ...
                    'callback',@(u,evnt)chgResDisp(E,k,'R'));
                s.Y=uimenu(m,'label','Y', ...
                    'callback',@(u,evnt)chgResDisp(E,k,'Y'));
                s.V1=uimenu(m,'label','V (avg. all points)','separator','on', ...
                    'callback',@(u,evnt)chgResDisp(E,k,'V1'));
                s.V3=uimenu(m,'label','V (avg. middle points)', ...
                    'callback',@(u,evnt)chgResDisp(E,k,'V3'));
                s.G1=uimenu(m,'label','G (avg. all points)', ...
                    'callback',@(u,evnt)chgResDisp(E,k,'G1'));
                s.G3=uimenu(m,'label','G (avg. middle points)', ...
                    'callback',@(u,evnt)chgResDisp(E,k,'G3'));
                s.R1=uimenu(m,'label','R (avg. all points)', ...
                    'callback',@(u,evnt)chgResDisp(E,k,'R1'));
                s.R3=uimenu(m,'label','R (avg. middle points)', ...
                    'callback',@(u,evnt)chgResDisp(E,k,'R3'));
                s.Y1=uimenu(m,'label','Y (avg. all points)', ...
                    'callback',@(u,evnt)chgResDisp(E,k,'Y1'));
                s.Y3=uimenu(m,'label','Y (avg. middle points)', ...
                    'callback',@(u,evnt)chgResDisp(E,k,'Y3'));
                set(s.(E.resdisp{k}),'checked','on')
                E.menuitems.resdisp(k) = s;
            end
        end
    end
    
    % Anchors and vessels
    methods
        function i = anchornew(E,p)
            i = length(E.anchors)+1;
            anchor = struct('active',true, ...
                'p',p, ...
                'vslidx',[], ...
                'hl',anchordisplay(E,p,i) ...
                );
            if i==1, E.anchors=anchor; else E.anchors(i)=anchor; end
            if nargout==0, clear i, end
        end
        function i = vslnew(E,edge,extidx)
            i = length(E.vsls)+1;
            vsl = struct('active',true, ...
                'edge',edge, ... % edge can be multiple (column vector)
                'extidx',extidx, ...
                'hl',vsldisplay(E,edge(1),i), ...
                'result',blankresult(E) ...
                );
            if i==1, E.vsls=vsl; else E.vsls(i)=vsl; end
            E.anchors(extidx(1)).vslidx(end+1) = i;
            E.anchors(extidx(2)).vslidx(end+1) = i;
            if nargout==0, clear i, end
        end
        function anchordelete(E,i)
            if ~isempty(E.anchors(i).vslidx), error programming, end
            delete(E.anchors(i).hl)
            E.anchors(i) = struct('active',false, ...
                'p',[], ...
                'vslidx',[], ...
                'hl',[] ...
                );
        end
        function vsldelete(E,i)
            % inform user about which vessel is deleted (give the ID)
            fprintf('removing vessel %s\n',E.vsls(i).edge(1).id_edge)
            
            % delete
            delete(E.vsls(i).hl)
            extidx = E.vsls(i).extidx;
            E.vsls(i) = struct('active',false, ...
                'edge',[], ...
                'extidx',[], ...
                'hl',[], ...
                'result', blankresult(E) ...
                );
            
            % modify extremities
            for k=extidx(:)'
                E.anchors(k).vslidx = setdiff(E.anchors(k).vslidx,i);
                if isempty(E.anchors(k).vslidx), anchordelete(E,k), end
            end
            
            % clear lines if vessel was current vessel
            if isequal(i,E.curvsl)
                E.curvsl = [];
                cla(E.hb(1))
                cla(E.hb(2))
            end
            
            % clear local parameters
            set(E.localparinput,'string','')
        end
        function hl = anchordisplay(E,p,i)
            hl(1) = line(p(1),p(2),'parent',E.ha(1));
            hl(2) = line(p(1),p(2),'parent',E.ha(2));
            set(hl,'marker','.','markersize',18,'color','r', ...
                'userdata',i, ...
                'hittest','on','buttondownfcn',@(hp,evnt)anchorclick(E,i))
        end
        function hl = vsldisplay(E,edge,i)
            if ~isscalar(edge), error('edge must be singleton'), end
            
            p = edge.snake(:,ceil(length(edge.snake)/3));
            hl = zeros(2,4);
            for k=1:2
                hl(k,1)=line(edge.snake(1,:),edge.snake(2,:), ...
                    'parent',E.ha(k),'hittest','off');
                hl(k,2)=line(edge.tubea(1,:),edge.tubea(2,:),'linestyle','--', ...
                    'parent',E.ha(k),'hittest','off');
                hl(k,3)=line(edge.tubeb(1,:),edge.tubeb(2,:),'linestyle','--', ...
                    'parent',E.ha(k),'hittest','off');
                hl(k,4)=line(p(1),p(2),'parent',E.ha(k),'marker','o','markersize',13);
                hl(k,5)=text(double(edge.snake(1,1)),double(edge.snake(2,1)),num2str(i), ...
                    'parent',E.ha(k),'hittest','off');
            end
            vslcolor(E,hl,edge)
            set(hl(2,:),'col','k')
            set(hl(:,4),'userdata',i, ...
                'hittest','on','buttondownfcn',@(hp,evnt)vslclick(E,i), ...
                'uicontextmenu',E.vslmenu)
        end
        function vslcolor(E,hl,edge) %#ok<MANU>
            if ~isscalar(edge), error('edge must be singleton'), end
            [col col0] = getcolor(edge);
            set(hl(1,4),'col',col0)
            set(hl(1,[1:3 5]),'col',col)
        end
        function vslflag(E,i,str)
            
            e = E.vsls(i).edge; % handle object
            [e.flag] = deal(str);
            vslcolor(E,E.vsls(i).hl,e(1))
        end
        function toggleactive(E,i)
            
            e = E.vsls(i).edge; % handle object
            [e.active] = deal(~e(1).active);
            vslcolor(E,E.vsls(i).hl,e(1))
        end
        function vslbad(E,i)
            
            e = E.vsls(i).edge; % handle object
            [e.flag] = deal('bad');
            [e.active] = deal(false);
            vslcolor(E,E.vsls(i).hl,e(1))
        end
        function anchorclick(E,i)
            
            switch get(gcf,'SelectionType')
                case 'normal'
                    switch E.selpolyflag
                        case ''
                            E.selpolyflag = 'start';
                        case 'current'
                            E.selpolyflag = 'end';
                    end
                    selectpoly(E,i)
            end
        end
        function vslclick(E,i)
            switch get(gcf,'SelectionType')
                case 'normal'
                    linesdisplay(E,i)
            end
        end
    end
    
    % Axes callback
    methods
        function axesclick(E,ha)
            
            if ismember(ha,E.ha)
                switch get(gcf,'SelectionType')
                    case 'open'             % set new anchor
                        p = get(ha,'currentpoint'); p = p(1,1:2)';
                        i = anchornew(E,p);
                        if strcmp(E.selpolyflag,'current')
                            % end of edge
                            E.selpolyflag = 'end';
                            selectpoly(E,i)
                        else
                            % start an edge
                            E.selpolyflag = 'start';
                            selectpoly(E,i)
                        end
                    case 'normal'           % selecting an edge
                        if strcmp(E.selpolyflag,'current')
                            selectpoly(E)
                        end
                    case {'alt','extend'}   % cancel current action
                        if strcmp(E.selpolyflag,'current')
                            E.selpolyflag = 'cancel';
                            selectpoly(E)
                        end
                end
            elseif ismember(ha,E.hb)
                switch get(gcf,'SelectionType')
                    case 'open'             % follow a trajectory
                        % only applies to image display (see function
                        % linesdisplay)
                        k = find(E.hb==ha);
                        if length(E.resdisp{k})~=1, return, end
                        % current vessel
                        i = E.curvsl;
                        if isempty(i) || isempty(E.vsls(i).result.V), return, end
                        % go!
                        p = get(ha,'currentpoint'); p = p(1,1:2)';
                        for i_res=1:4
                            dum = reshape('VyGgRbYr',2,4);
                            f = dum(1,i_res);
                            col = dum(2,i_res);
                            V = E.vsls(i).result.(f);
                            if isempty(V), continue, end
                            [nx nt] = size(V);
                            t = round(p(1)); x = p(2);
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
                            xdata = [xa(end:-1:1) x0 xb];
                            tdata = t0 + (-length(xa):length(xb));
                            hl = [];
                            for k2=1:5
                                if length(E.resdisp{k2})==1
                                    hl(end+1) = line(tdata,xdata, ...
                                        'color',col, ...
                                        'parent',E.hb(k2), ...
                                        'tag','traj'); %#ok<AGROW>
                                end
                            end
                            set(hl,'tag','trajectory','hittest','on', ...
                                'buttondownfcn',@(h,evnt)trydelete(hl,E.hf))
                        end
                end
            end
        end
    end
    
    % Edge segmentation
    methods
        function selectpoly(E,i)
            switch E.selpolyflag
                case 'start'
                    p = E.anchors(i).p;
                    E.tmp.poly = p;
                    E.tmp.hl = initpoly(E);
                    E.tmp.extidx = i;
                    set(gcf,'windowbuttonmotionfcn',@(h,evnt)drawpoly(E))
                    E.selpolyflag = 'current';
                case 'current'
                    p = get(gca,'currentpoint'); p = p(1,1:2)';
                    E.tmp.poly(:,end+1) = p;
                case 'end'
                    % finish current polygon selection
                    set(gcf,'windowbuttonmotionfcn','')
                    E.selpolyflag = '';
                    delete(E.tmp.hl)
                    % set new vessel
                    p = E.anchors(i).p;
                    E.tmp.poly(:,end+1) = p;
                    E.tmp.extidx(2) = i;
                    setvessel(E)
                case 'cancel'
                    % cancel current polygon selection
                    set(gcf,'windowbuttonmotionfcn','')
                    delete(E.tmp.hl)
                    E.selpolyflag = '';
                    % delete the starting anchor if it is not connected
                    i = E.tmp.extidx;
                    if isempty(E.anchors(i).vslidx), anchordelete(E,i), end
                otherwise
                    error programming
            end
        end
        function hl = initpoly(E)            
            poly = E.tmp.poly;
            hl = [0 0];
            for i=1:2
                hl(i) = line(poly(1,:),poly(2,:),'parent',E.ha(i), ...
                    'color',[.5 .5 .5]);
            end
        end
        function drawpoly(E)            
            p = get(gca,'currentpoint'); p = p(1,1:2)';
            poly = [E.tmp.poly p];
            set(E.tmp.hl,'xdata',poly(1,:),'ydata',poly(2,:))
        end
        function setvessel(E)           
            poly = E.tmp.poly;
            
            % automatic segmentation
            edge = fast_segment(E.S.mask,poly,E.parameters.segment,E.ha);
            edge.flag = '';
            
            % duplicate edge if multiple data
            global Y
            if iscell(Y)
                for k=2:numel(Y), edge(k,1) = copy(edge(1)); end
            end
            
            % new vessel
            i = vslnew(E,edge,E.tmp.extidx);
            
            % display lines
            linesdisplay(E,i)
        end
        function upsidedown(E,i)
            vsl = E.vsls(i);
            delete(E.vsls(i).hl)
            E.vsls(i).hl = [];
            
            % inversion
            edge = vsl.edge;
            edge = fast_edge(edge.snake(:,end:-1:1),edge.halfd(:,end:-1:1),edge.dx);
            edge.flag = '';
            
            % duplicate edge if multiple data
            global Y
            if iscell(Y)
                for k=2:numel(Y), edge(k,1) = copy(edge(1)); end
            end
            
            % change vessel
            E.vsls(i).edge = edge;
            E.vsls(i).hl = vsldisplay(E,edge(1),i);
            E.vsls(i).result = blankresult(E); % force re-compute
            
            % display lines
            linesdisplay(E,i)
        end
        function resetvessel(E,i)            
            vsl = E.vsls(i);
            poly = [E.anchors(vsl.extidx(1)).p vsl.edge(1).snake ...
                E.anchors(vsl.extidx(2)).p];
            delete(E.vsls(i).hl)
            E.vsls(i).hl = [];
            
            % automatic segmentation
            edge = fast_segment(E.S.mask,poly,E.parameters.segment,E.ha);
            edge.flag = '';
            
            % duplicate edge if multiple data
            global Y
            if iscell(Y)
                for k=2:numel(Y), edge(k,1) = copy(edge(1)); end
            end
            
            % change vessel
            E.vsls(i).edge = edge;
            E.vsls(i).hl = vsldisplay(E,edge(1),i);
            E.vsls(i).result = blankresult(E); % force re-compute
            
            % display lines
            linesdisplay(E,i)
        end
        function refilter(E)
            i = E.curvsl;
            if isempty(i), return, end
            E.vsls(i).result.I = [];
            linesdisplay(E,i)
        end
    end
    
    % Lines image
    methods
        function linesdisplay(E,i)            
            % mark current edge (and un-mark previous)
            if ~isempty(E.curvsl)
                set(E.vsls(E.curvsl).hl,'linewidth',1)
            end
            E.curvsl = i;
            set(E.vsls(i).hl,'linewidth',2)
            drawnow
            
            % display note
            set(E.noteinput,'string',E.vsls(i).edge(1).note)
            
            % update local parameters
            set(E.localparinput,'string', ...
                fn_struct2str(E.vsls(i).edge(E.kY).estpar_local))
            
            % lines data
            global Y
            if isempty(E.vsls(i).result.I0)
                % interpolation
                edge = E.vsls(i).edge(E.kY);
                if iscell(Y)
                    y = fast_interptube(Y{E.kY},edge);
                else
                    y = fast_interptube(Y,edge);
                end
                % store
                E.vsls(i).result.I0 = y;
                % fake true speed
                if status(E.S,'fake')
                    x = E.S.faketruespeed{E.kY};
                    snake = E.vsls(i).edge(E.kY).snake(1,:);
                    E.vsls(i).result.Y = interp1(1:E.S.nx,x,snake);
                end
            end
            
            % filtering
            if isempty(E.vsls(i).result.I)
                % filter
                y = E.vsls(i).result.I0;
                str = E.parameters.linefilter;
                for k=1:length(str)
                    eval(str{k})
                end
                % store
                E.vsls(i).result.I = y;
            end
            
            % do not compute results if they don't exist
            
            % memorize the current axis
            ax = [];
            for k=1:5
                if fn_ismemberstr(E.resdisp{k},{'I','J','V','G'})
                    ax = axis(E.hb(1));
                    if ax(3)==0, ax=[]; end % this happens at 'fast_edgesmouse' start
                    break
                end
            end
            
            % display parameters
            par = fn_str2struct(E.parameters.resdisp);

            % display lines and results
            for k=1:5
                f = E.resdisp{k}; % f is a letter + possibly a number
                y = E.vsls(i).result.(f(1)) * par.fact;
                if ~any(f(1)=='IJ') && mean(y(:))<0, y=-y; end % make positive speed negative
                if isempty(y), y = zeros(size(E.vsls(i).result.I0)); end
                if length(f)==1
                    % IMAGE DISPLAY
                    [np nt] = size(y);
                    ycut = y(ceil(np/3):ceil(2*np/3),ceil(nt/10):ceil(9*nt/10));
                    clip = [min(ycut(:)) max(ycut(:))];
                    if ~diff(clip), clip = clip + [-1 1]; end
                    imagesc(y,'parent',E.hb(k),clip)
                    if ~isempty(ax), fn_imvalue('chgx',ax(1:2),E.hb(k)), end
                    set(E.hb(k),'ydir','normal')
                else
                    % PLOT DISPLAY
                    % the number indicates how to average in the spatial
                    % dimention
                    nx = size(y,1);
                    switch f(2)
                        case '1'
                            idx = 1:nx;
                        case '3'
                            idx = round(nx/3+1):round(2*nx/3);
                        otherwise
                              error programming
                    end
                    % time courses
                    y = abs(mean(y(idx,:),1));
                    if par.L
                        y = filtx(y,par.L,'lm');
                    end
                    plot(y,'parent',E.hb(k))
                    axis(E.hb(k),'tight')
                end
            end
            
            % display average estimated speed
            V = E.vsls(i).result.V;
            if isempty(V)
                str = '';
            else
                G = E.vsls(i).result.G;
                if isempty(G)
                    str = sprintf('mean speed: %.1f',mean(V(:))*par.fact);
                else
                    str = sprintf('mean speed: %.1f (Gabor: %.1f)',mean(V(:))*par.fact,mean(G(:))*par.fact);
                end
            end
            set(E.statusbar,'string',str)
            
        end
        function linestrack(E)
            % get lines data
            i = E.curvsl;
            if isempty(i), return, end
            result = E.vsls(i).result;
            if isempty(result.I0) || isempty(result.I), return, end
            y = result.I;
            
            % estimation parameters
            par = fn_str2struct(E.parameters.estpar);
            par = fn_structmerge(par,E.vsls(i).edge(E.kY).estpar_local);
            
            % compute
            set(findobj(E.hb(1),'tag','traj'),'color','r')
            set(E.hf,'pointer','watch')
            fn_progress('in',E.statusbar)
            res = fast_estimateflow(y,par,E.hb(2:3));
            fn_progress screen
            E.vsls(i).result = fn_structmerge(E.vsls(i).result,res);
            set(E.hf,'pointer','arrow')
            
            % display
            linesdisplay(E,i)
        end
        function chgResDisp(E,k,resname)
            % update marks
            set(E.menuitems.resdisp(k).(E.resdisp{k}),'checked','off')
            set(E.menuitems.resdisp(k).(resname),'checked','on')
            % update value
            E.resdisp{k} = resname;   
            % update 'lines' display
            i = E.curvsl;
            if isempty(i), return, end
            linesdisplay(E,i)
        end 
        function chgresdisppar(E,str)
            % update options
            E.parameters.resdisp = str;
            % update 'lines' display
            i = E.curvsl;
            if isempty(i), return, end
            linesdisplay(E,i)
        end
        function s = blankresult(E) %#ok<MANU>
            s = struct('I0',[],'I',[],'J',[],'V',[],'G',[],'R',[],'Y',[]);
        end
    end
    
    % Edit parametrs
    methods
        function chgFakeIdx(E,dim,idx)
            global Y
            if ~iscell(Y), error programming, end
            % set the appropriate kY
            s = size(Y);
            ijk = fn_indices(s,E.kY);
            set(E.menuitems.fake{dim}(ijk(dim)),'checked','off')
            ijk(dim) = idx;
            set(E.menuitems.fake{dim}(ijk(dim)),'checked','on')
            E.kY = fn_indices(s,ijk);
            % re-display lines
            i = E.curvsl;
            if isempty(i), return, end
            E.vsls(i).result = blankresult(E); % force re-compute
            linesdisplay(E,i)
        end
        function chgFakeData(E,doallidx)
            global Y
            % do all indices
            if doallidx
                kys = 1:numel(Y);
            else
                kys = E.kY;
            end
            % recompute
            fn_progress('create fake data',length(kys))
            s = size(Y); s(end+1:4)=1;
            for i=1:length(kys)
                ky = kys(i);
                fn_progress(i)
                ijk = fn_indices(s,ky);
                createfake(E.S,ijk(1),ijk(2),ijk(3),ijk(4))
            end
            % re-display lines
            i = E.curvsl;
            if isempty(i), return, end
            E.vsls(i).result = blankresult(E); % force re-compute
            linesdisplay(E,i)
        end
        function chgsegmentpar(E,str)
            E.parameters.segment = cellstr(str);
            if ~isempty(E.curvsl) && E.options.doresponsiveest
                resetvessel(E,E.curvsl)
            end
        end
        function chglinefilter(E,str)
            E.parameters.linefilter = cellstr(str);
            if ~isempty(E.curvsl) && E.options.doresponsiveest
                refilter(E)
            end
        end
        function chglinetrackglobal(E,str)       
            E.parameters.estpar = cellstr(str);
            if ~isempty(E.curvsl) && E.options.doresponsiveest
                linestrack(E)
            end
        end
        function chglinetracklocal(E,str) 
            if isempty(E.curvsl)
                set(E.localparinput,'string','')
            else
                par = fn_str2struct(str);
                if isfield(par,'method')
                    errordlg('''method'' cannot be overwritten')
                    return
                end
                E.vsls(E.curvsl).edge(E.kY).estpar_local = par;
                if E.options.doresponsiveest
                    linestrack(E)
                end
            end
        end
    end
    
    % Save edges
    methods
        function savebutton(E)
            e = [E.vsls.edge]; % note that edge is a handle object
            
            % save the estimation parameters into e
            [e.filterpar] = deal(E.parameters.linefilter);
            par = fn_str2struct(E.parameters.estpar);
            for i=1:numel(e)
                e(i).estpar = fn_structmerge(par,e(i).estpar_local);
            end
            
            % save
            E.S.edges = e;
            E.S.parameters = fn_structmerge(E.S.parameters,E.parameters);
            E.S.id_res = [];
            saveproject(E.S)
        end
    end
    
    % Misc
    methods
        function access(E)
            keyboard
        end
    end
end


%---
function autoclip(ha,x)

ax = axis(ha);
ii = round(ax(1)):round(ax(2)-1e-6);
jj = round(ax(3)):round(ax(4)-1e-6);
x = x(ii,jj);
clip = [min(x(:)) max(x(:))];
if diff(clip)==0, clip = clip+[-1 1]; end
set(ha,'clim',clip)
end

%---
function trydelete(hl,hf)
switch get(hf,'selectionType')
    case 'normal'
        delete(hl(ishandle(hl)))
    otherwise
        delete(findobj(hf,'tag','trajectory'))
end
end
        

%---
function copyDisplay(ha)
% copy display in new figure

str = inputdlg('name?');
hf = figure('numbertitle','off','name',str{1});
c = copyobj(ha,hf);
set(c,'pos',[.13 .11 .775 .815])

end




