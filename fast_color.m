classdef fast_color < interface
    
    properties (SetAccess='private')
        S
        ntime   % S.nt if average over trials, S.nexp if average over time
        a                   % sparse array to transform vessel result into movie
        indices             % indices in image where there are vessels
        x = cell(1,2);      % vessel results in separate cells
        xall = cell(1,2);   % vessel results alltogether
        CS      % reference image
        CSidx 	% reference image values at indices
        Xcol    % display controls
        Xslid   % time slider control
        Xplay   % play speed control
        fps
    end
    properties (Dependent, SetAccess='private')
        colpar
        kframe
    end
    properties (Access='private')
        im
        hline
        hnumbers
    end
    
    % INIT
    methods
        function C = fast_color(S)
            if nargin==0, S = evalin('base','S'); end
            
            % initialization of the parent interface
            defaultoptions = struct;
            % parameters for the 2 displays
            defaultoptions.parcol = struct( ...
                'data',     {'volume','result'}, ...
                'algo',     'track', ...
                'trials',   'all', ...
                'res',      'stim/rest', ...
                'fzsub',    [], ...
                'tsmooth',  [], ...
                'xsmooth',  [],  ...
                'clip',     [-1 1], ...
                'cmap',     'hot', ...
                'overlay',  .5 ...
                );
            % speed control
            defaultoptions.parplay = struct('playspeed', 5, ...
                'framestep', 1);
            C = C@interface(287,'VESSEL COLOR',defaultoptions);

            % checks
            if isempty(S.edges), errordlg('no edges'), return, end
            if isempty(S.nexp), errordlg('no project to display'), return, end
                        
            % some initializations
            set(0,'defaultfigurecolormap',gray(256))
            
            % project info
            C.S = S;
            C.ntime = fn_switch(S.parameters.resultavgtrials,S.nt,S.nexp);
            
            % graphic objects
            set(C.hf,'defaultuicontrolfontsize',10)
            C.grob.hp(1:2) = [uipanel uipanel];
            C.grob.haim(1:2) = [axes axes];
            C.grob.hagr(1:2) = [axes axes];
            C.grob.hplay = uipanel;
            C.grob.hslid = uipanel;
            
            % end of interface initialization
            % calls init_menus -> create a menu in the case of 'fake' data
            interface_end(C), drawnow
            S.closefig = union(S.closefig,C.hf);
            
            % indices operations
            C.a = v2iTransform(C);
            C.indices = find(any(C.a,2));
            C.a = C.a(C.indices,:);
            
            % init graphic display
            init_display(C);
            
            % init display options
            init_colpar(C)
            
            % init slider and control of speed
            init_player(C)
            
            % load results
            loadresult(C,1)
            loadresult(C,2)
            
            % go
            displaygraph(C)
            movieshow(C)

        end
        function init_menus(C)
            init_menus@interface(C)
            uimenu(C.menus.interface, ...
                'label','make current parameters default','separator','on', ...
                'callback',@(u,evnt)makedefaultpar)
            function makedefaultpar
                C.options.parcol = C.colpar;
                C.options.parplay = C.Xplay.s;
                saveoptions(C)
            end
        end
        function init_colpar(C)
            % specification for display parameters
            colparspec = struct( ...
                'data',     {{'volume','result','width'}}, ...
                'algo',     {{'top','track','radon'}}, ...
                'trials',   {{'all','odd','even'}}, ...
                'res',      {{'all','stim/rest','stim','rest'}}, ...
                'fzsub',    'xslider 0 200 1', ...
                'tsmooth',  'xslider 0 50 1', ...
                'xsmooth',  'xslider 0 20 1',  ...
                'clip',     'clip', ...
                'cmap',     {{'gray','jet','hot','mapclip','mapcliphigh','signcheck','green'}}, ...
                'overlay',  'slider 0 1' ...
                );
            % parameters controls
            C.Xcol    = fn_control(C.options.parcol(1),@(s)updcolpar(C,1), ...
                colparspec,C.grob.hp(1));
            C.Xcol(2) = fn_control(C.options.parcol(2),@(s)updcolpar(C,2), ...
                colparspec,C.grob.hp(2));
        end
        function init_display(C)
            % frame display and vessel numbers
            m = min(C.S.CS(:));
            M = max(C.S.CS(:));
            C.CS = (C.S.CS-m)/(M-m);
            C.CS = C.CS';
            C.CS = repmat(C.CS(:),[1 3]);
            C.CSidx = C.CS(C.indices,:);
            edges = C.S.edges([C.S.edges.active]);
            ne = length(edges);
            C.hnumbers = zeros(ne,2);
            for k=1:2
                ha = C.grob.haim(k);
                % frame
                fr = reshape(C.CS,[C.S.ny C.S.nx 3]);
                C.im(k) = image(fr,'parent',ha);
                set(ha,'xtick',[],'ytick',[])
                % vessel numbers
                for i=1:length(edges)
                    C.hnumbers(i,k) = text( ...
                        double(edges(i).snake(1,1)),double(edges(i).snake(2,1)),num2str(i), ...
                        'parent',ha);
                end
            end
            % graph - nothing to do
        end
        function init_player(C)
            % time slider
            C.Xslid = fn_slider(C.grob.hslid,'mode','point', ...
                'minmax',[1 C.ntime],'value',1, ...
                'inc',1/(C.ntime-1),'width',1/(C.ntime-1)', ...
                'callback',@slide);
            function slide(u,evnt) %#ok<INUSD>
                displayframe(C)
                set(C.hline,'xdata',C.kframe([1 1]))
            end
            % speed control
            spec = struct('playspeed', 'xlogslider 0 2 %.0f', ...
                'framestep', 'slider 1 50 1');
            C.fps = C.options.parplay.playspeed;
            C.Xplay = fn_control(C.options.parplay, ...
                @(u,evnt)updplay(C),spec,C.grob.hplay);
        end
        function movieshow(C)
            while ishandle(C.hf) && ~isempty(C.fps)
                tic
                % automatic display update
                C.kframe = 1+mod(C.kframe+C.Xplay.framestep-1,C.ntime); 
                if ishandle(C.hf) && ~isempty(C.fps), pause(1/C.fps-toc), end
            end
        end
    end
    
    % GET/SET
    methods
        function k = get.kframe(C)
            k = C.Xslid.value;
        end
        function set.kframe(C,k)
            C.Xslid.value = k;
            displayframe(C)
            if ~ishandle(C.hf), return, end
            set(C.hline,'xdata',C.kframe([1 1]))
        end
        function s = get.colpar(C)
            for k=1:2
                s(k) = C.Xcol(k).s; %#ok<AGROW>
            end
        end
    end
    
    % vessel display tools
    methods (Access='private')
        function updcolpar(C,k)
            % need to reload results?
            if any(fn_ismemberstr(C.Xcol(k).changedfields, ...
                    {'data','algo','trials','res','fzsub','tsmooth','xsmooth'}))
                % operation can be long
                set(C.hf,'pointer','watch'), drawnow
                loadresult(C,k)
                set(C.hf,'pointer','arrow')
                displaygraph(C,k)
            end
            % display update
            displayframe(C,k)
        end
        function updplay(C)
            oldfps = C.fps;
            C.fps = C.Xplay.playspeed;
            if isempty(oldfps) && ~isempty(C.fps), movieshow(C), end
        end
        function loadresult(C,k)
            % result parameters
            par = C.colpar(k);
            resdisp = struct( ...
                'data',     par.data, ...
                'algo',     par.algo, ...
                'trial',    'average', ...
                'res',      par.res, ...
                'oddeven',  par.trials, ...
                'points',   -par.xsmooth, ...
                'smoothing',par.tsmooth, ...
                'fzsub',    par.fzsub ...
                );
            if ~C.S.parameters.resultavgtrials, resdisp.trial='timeavg'; end
            if isempty(par.xsmooth), resdisp.points='all'; end
            if isempty(par.fzsub),   resdisp.fzsub=0; end
            if strcmp(resdisp.data,'width'), resdisp.data = 'sigma'; end
            % get result
            C.x{k} = getresult(C.S,resdisp);
            % handle vessels without result (not computed yet)
            kedges = find([C.S.edges.active]);
            ne = length(kedges);
            if ne~=length(C.x{k}), error programming, end
            for i=1:ne
                if isempty(C.x{k}{i})
                    if isempty(par.xsmooth), np=1; else np=C.S.edges(kedges(i)).np; end
                    C.x{k}{i} = zeros(np,C.ntime);
                end
            end
            % average on all points: repmat
            if isempty(par.xsmooth)
                for i=1:length(C.x{k})
                    y = C.x{k}{i};
                    C.x{k}{i} = repmat(y(:)',C.S.edges(kedges(i)).np,1);
                end
            end
            % concatenate all vessels together
            C.xall{k} = double(cat(1,C.x{k}{:}));
        end
        function displayframe(C,klist)
            if nargin<2, klist=1:2; end
            for k=klist
                if ~ishandle(C.hf), return, end
                % parameters
                p = C.Xcol(k).s;
                cmap = feval(p.cmap,256);
                % movie frame
                y = C.a * C.xall{k}(:,round(C.kframe));
                y = (y-p.clip(1)) * (256/(p.clip(2)-p.clip(1)));
                y = min(256,max(1,ceil(y)));
                y = cmap(y(:),:);
                fr = (1-p.overlay)+p.overlay*C.CS;
                fr(C.indices,:) = y;
                fr = reshape(fr,[C.S.ny C.S.nx 3]);
                set(C.im(k),'cdata',fr);
                drawnow
            end
        end
        function displaygraph(C,klist)
            if nargin<2, klist=1:2; end
            for k=klist
                if ~ishandle(C.hf), return, end
                % average over vessels
                ne = length(C.x{k});
                data = zeros(C.ntime,ne);
                for i=1:ne, data(:,i)=mean(C.x{k}{i},1); end
                % display
                ha = C.grob.hagr(k);
                fn_eegplot(data,'parent',ha,'5STD','num')
                axis tight
                ax = axis(ha);
                C.hline(k) = line(C.kframe([1 1]),ax(3:4),'parent',ha, ...
                    'color','k');
                drawnow
            end
            
        end
        function a = v2iTransform(C,e)
            % get the sparse array to interpolate vessel data into an image
            S = C.S;
            % all edges
            if nargin<2
                edges  = S.edges([S.edges.active]);
                np = sum([edges.np]);
                a = sparse(S.nx*S.ny,0);
                for i=1:numel(edges)
                    a = [a v2iTransform(C,edges(i))]; %#ok<AGROW>
                end
                return
            end
            % operation
            x=double([e.tubea e.tubeb(:,end:-1:1)]);
            mask=poly2mask(x(2,:),x(1,:),S.nx,S.ny);
            [ii jj] = find(mask); % column vectors
            dd = fn_add(ii,-e.snake(1,:)).^2 + fn_add(jj,-e.snake(2,:)).^2;
            [dum kk] = min(dd,[],2); % indices in vessels corresponding to minimal distance
            %a = sparse(ii+S.nx*(jj-1),kk,1,S.nx*S.ny,e.np);
            a = sparse(jj+S.ny*(ii-1),kk,1,S.nx*S.ny,e.np);
        end
    end
    
    % Misc
    methods
        function access(C)
            keyboard
        end
    end
end

%---
function trydelete(hl)

delete(hl(ishandle(hl)))

end
    
%---
