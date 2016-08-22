classdef fast_displayresult < interface
    
    properties (SetAccess='private')
        S
        hl
        e
        volumeroi
        dx
    end
    properties
        dispvolroi = false;
        chgtrials = false;
        dxcorr = true;
    end
    properties (Dependent)
        kedge
        ktrial
    end
    properties (Access='private')
        kedg
        ktri
        kY
        trimg
        menuitems
        resdispoptions
        inputnote
    end
    
    methods
        function R = fast_displayresult(S,nexpmax)
            if nargin==0, S = evalin('base','S'); end
            if nargin<2, nexpmax = S.nexp; end
            
            % initialization of the parent interface
            defaultoptions.resdisp = struct( ...
                'data',     {'volume','volume','flow','flow'}, ...
                'trials',   {'odd+even','all','odd+even','all'}, ...
                'points',   {'middle','middle','middle','middle'}, ...
                'res',      {'stim+rest','stim/rest','stim+rest','stim/rest'}, ...
                'fzsub',    {0,5,0,5}, ...
                'smoothing',{0,30,0,30} ...
                );
            R = R@interface(S.hf+2,'RESULT',defaultoptions);
%             for i=1:4
%                 if ~any(strcmp(R.options.resdisp(i).data,[S.parameters.resultnames {'compare'}]))
%                     R.options.resdisp(i).data = 'volume';
%                 end
%             end
                    
            % checks
            if isempty(S.edges), errordlg('no edges'), return, end
            if isempty(S.nexp), errordlg('no project to display'), return, end
                        
            % some initializations
            fn_imvalue image
            set(0,'defaultfigurecolormap',gray(256))
            
            % project info
            R.S = S;
            R.ktri = 1;
            R.kY = 1;
            
            % vasculature image
            R.grob.haim = axes;
            display_image(R)
            
            % check that all vessels where interpolated at the same spacing
            R.dx = [R.S.edges.dx];
            if any(diff(R.dx)), error('edges where not interpolated with the same spacing'), end
            R.dx = R.dx(1);
            
            % parameters of the coregistration
            R.grob.haiso = axes;
            display_iso(R)
            
            % result display
            R.grob.hares = [axes axes axes axes];
            
            % estimation result on single trial
            R.grob.halin = [axes axes axes];
            R.grob.hulinshift = uicontrol('style','popupmenu', ...
                'string',{'no shift', ...
                '1pix','2pix','3pix','4pix','5pix','6pix','7pix','8pix','9pix','10pix', ...
                '-1pix','-2pix','-3pix','-4pix','-5pix','-6pix','-7pix','-8pix','-9pix','-10pix'}, ...
                'value',1, ...
                'callback',@(u,evt)display_trialest(R));

            % buttons
            R.grob.hu(1) = uicontrol('string','PRECOMP', ...
                'callback',@(hu,evnt)precomp(R.S));
                        
            R.grob.hu(2) = uicontrol('string','MOVIE', ...
                'callback',@movie);
            if status(S,'fake'), set( R.grob.hu(2),'visible','off'), end
            
            function movie(u,evnt) %#ok<INUSD>
                set(R.hf,'pointer','watch'), drawnow
                Z = gettrial(R.S,R.ktrial,'int16');
                set(R.hf,'pointer','arrow'), drawnow
                fn_movie(Z)
            end
           
            % note and info
            R.grob.hu(3) = uicontrol('style','text','fontsize',8, ...
                'string','note:');
            R.grob.hu(4) = uicontrol('style','edit','max',1, ...
                'fontsize',8,'horizontalalignment','left', ...
                'backgroundcolor','w', ...
                'callback',@(hu,evnt)set(R.e,'note',get(hu,'string')));
            R.inputnote = R.grob.hu(4);
            R.grob.huinfo = uicontrol('style','text','fontsize',8, ...
                'horizontalalignment','left');
            
            % click in window
            set(R.grob.hf,'windowbuttondownfcn',@(hf,evnt)axesclick(R,gca))
            
            % that's it for the graphic initializations, position windows
            interface_end(R)

%             % automatic displays upon definition of kedge!
%             R.kedge = 1;
        end
        function init_menus(R)
            init_menus@interface(R)
            if isfield(R.menus,'fake'), delete(R.menus.fake), end
            if isfield(R.menus,'resdisp'), delete(R.menus.resdisp), end
            if isfield(R.menus,'volumeroi'), delete(R.menus.volumeroi), end
            S = R.S;
            
            % correct estimation results by interpolating distance?
            R.menuitems.dxcorr = uimenu(R.menus.interface, ...
                'label','correct results by spacing','separator','on', ...
                'checked',fn_switch(R.dxcorr), ...
                'callback',@(u,evnt)set(R,'dxcorr',~R.dxcorr));
            % allow changing the trial selection?
            R.menuitems.dxcorr = uimenu(R.menus.interface, ...
                'label','change trial selection', ...
                'checked',fn_switch(R.chgtrials), ...
                'callback',@(u,evnt)fn_setpropertyandmark(R,'chgtrials',u,'switch'));
            
            % fake data: create a menu to choose which data to use
            if status(S,'fake')
                s = S.fake;
                m = uimenu(R.hf,'label','fake parameters');
                items = cell(1,4);
                for i=1:s.nst
                    mi=uimenu(m,'label',['step ' num2str(s.step(i))], ...
                        'callback',@(u,evnt)chgY(R,1,i));
                    if i==1, set(mi,'checked','on'), end
                    items{1}(end+1)=mi;
                end
                for i=1:s.nth
                    mi=uimenu(m,'label',['thresh ' num2str(s.thresh(i))], ...
                        'callback',@(u,evnt)chgY(R,2,i));
                    if i==1, set(mi,'checked','on','separator','on'), end
                    items{2}(end+1)=mi;
                end
                for i=1:s.nji
                    mi=uimenu(m,'label',['jitter ' num2str(s.jitter(i))], ...
                        'callback',@(u,evnt)chgY(R,3,i));
                    if i==1, set(mi,'checked','on','separator','on'), end
                    items{3}(end+1)=mi;
                end
                for i=1:s.nsn
                    mi=uimenu(m,'label',['shotnoise ' num2str(s.shotnoise(i))], ...
                        'callback',@(u,evnt)chgY(R,4,i));
                    if i==1, set(mi,'checked','on','separator','on'), end
                    items{4}(end+1)=mi;
                end
                R.menus.fake = m;
                R.menuitems.fake = items;
            end
            
            % switch between different estimation parameter sets
            
            % list of options for result display
%             R.resdispoptions = struct( ...
%                 'data',     {[{'volume'} R.S.parameters.resultnames {'compare'}]}, ...
%                 'trials',   {{'all','current','odd+even','odd','even'}}, ...
%                 'res',      {{'all','stim/rest','stim+rest','stim','rest'}}, ...
%                 'fzsub',    {{0,5,50,200}}, ...
%                 'smoothing',{{0,5,10,20,30}}, ...
%                 'points',   {{'all','middle','2/5','3/5','4/5'}} ...
%                 );
            % menus to select what to see in the global results graphs
            m = zeros(1,4);
            items = struct;
%             for k = 1:4
%                 m(k) = uimenu(R.hf,'label',['result display ' num2str(k)]);
%                 F = fieldnames(R.resdispoptions);
%                 for i=1:length(F)
%                     f = F{i};
%                     for j = 1:length(R.resdispoptions.(f))
%                         val = R.resdispoptions.(f){j};
%                         str = fn_switch(isnumeric(val),num2str(val),val);
%                         items(k).(f)(j) = uimenu(m(k),'label',[f ': ' str], ...
%                             'callback',@(u,evnt)chgResDisp(R,k,f,val));
%                         if i~=1 && j==1
%                             set(items(k).(f)(j),'separator','on')
%                         end
%                         if isequal(val,R.options.resdisp(k).(f))
%                             set(items(k).(f)(j),'checked','on')
%                         end
%                     end
%                 end
%             end
            R.menus.resdisp = m;
            R.menuitems.resdisp = items;
            
            % menu for superimposing volume response in a region of
            % interest
            m = uimenu(R.hf,'label','volume in ROI');
            items = struct;
            items.define = uimenu(m,'label','define ROI...', ...
                'callback',@(u,evnt)volumeROI(R));
            items.show = uimenu(m,'label','superimpose to results', ...
                'callback',@(u,evnt)set(R,'dispvolroi',~fn_switch(get(u,'checked'))));
            R.menus.volumeroi = m;
            R.menuitems.volumeroi = items;
        end
    end
    
    % GET/SET
    methods
        function i = get.kedge(R)
            i = R.kedg;
        end
        function k = get.ktrial(R)
            k = R.ktri;
        end
        function set.kedge(R,i)
            if i==R.kedg, return, end
            set(R.hf,'pointer','watch'), drawnow
            set(R.hl(R.kedg),'linewidth',1)
            set(R.hl(i),'linewidth',3)
            R.kedg = i;
            R.e = R.S.edges(R.kY,i);
            set(R.inputnote,'string',R.e.note), drawnow
            if isfield(R.e.estpar,'v0'), v0=R.e.estpar.v0; else v0=[]; end
            vmax = R.e.estpar.vmax; if isscalar(vmax), vmax=vmax*[-1 1]; end
            set(R.grob.huinfo,'string',sprintf('v0: %.1f\nvmax: %.1f %.1f',v0,vmax(1),vmax(2)))
            display_result(R), drawnow
            set(R.grob.hulinshift,'value',1) % cancel line shift
            display_trialest(R)            
            set(R.hf,'pointer','arrow')
        end
        function set.ktrial(R,k)
            if k==R.ktri, return, end
            R.ktri = k;
            % update display of current trial result
            display_trialest(R)
            % update display of global results if, in fact, they show only
            % the current trial
            displayscurrent = strcmp({R.options.resdisp.trials},'current');
            if any(displayscurrent)
                display_result(R,find(displayscurrent))
            end
        end
        function chgY(R,dim,idx)
            % set the appropriate kY
            s = [R.S.fake.nst R.S.fake.nth R.S.fake.nji R.S.fake.nsn];
            ijk = fn_indices(s,R.kY);
            if idx==ijk(dim), return, end
            set(R.menuitems.fake{dim}(ijk(dim)),'checked','off')
            ijk(dim) = idx;
            set(R.menuitems.fake{dim}(ijk(dim)),'checked','on')
            R.kY = fn_indices(s,ijk);
            % re-display results
            set(R.hf,'pointer','watch'), drawnow
            R.e = R.S.edges(R.kY,R.kedg);
            display_result(R)
            display_trialest(R)
            set(R.hf,'pointer','arrow')
        end
        function chgResDisp(R,k,f,val)
            % update check marks
            oldval = R.options.resdisp(k).(f);
            if isnumeric(oldval)
                j = find(oldval==cell2mat(R.resdispoptions.(f)));
            else
                j = strcmp(oldval,R.resdispoptions.(f));
            end
            set(R.menuitems.resdisp(k).(f)(j),'checked','off');
            if isnumeric(val)
                j = find(val==cell2mat(R.resdispoptions.(f)));
            else
                j = strcmp(val,R.resdispoptions.(f));
            end
            set(R.menuitems.resdisp(k).(f)(j),'checked','on');
            % update value
            R.options.resdisp(k).(f) = val;
            saveoptions(R)
            % update display
            display_result(R,k)
        end
        function set.dxcorr(R,b)
            if R.dxcorr==b, return, end
            R.dxcorr = b;
            set(R.menuitems.dxcorr,'checked',fn_switch(b)); %#ok<MCSUP>
            display_result(R);
        end
        function set.dispvolroi(R,b)
            if R.dispvolroi==b, return, end
            R.dispvolroi = b;
            % update menu mark
            set(R.menuitems.volumeroi.show,'checked',fn_switch(b))
            % set ROI if does not exist
            if isempty(R.volumeroi), volumeROI(R), end
            % update result display
            display_result(R)
        end
    end

    % display
    methods
        function display_image(R)
            imagesc(R.S.CS','parent',R.grob.haim)
            set(R.grob.haim,'xtick',[],'ytick',[])
            ne = size(R.S.edges,2);
            R.hl = zeros(1,ne);
            for i=1:ne
                edg = R.S.edges(1,i);
                R.hl(i) = line(edg.snake(1,:),edg.snake(2,:), ...
                    'parent',R.grob.haim,'col',getcolor(edg), ...
                    'hittest','on','busyaction','cancel', ...
                    'buttondownfcn',@(hl,evnt)set(R,'kedge',i));
                text(double(edg.snake(1,1)),double(edg.snake(2,1)),num2str(i), ...
                    'col',getcolor(edg),'hittest','off')
                line(edg.snake(1,1),edg.snake(2,1), ...
                    'marker','.', ...
                    'parent',R.grob.haim,'col',getcolor(edg), ...
                    'hittest','off', ...
                    'buttondownfcn',@(hl,evnt)set(R,'kedge',i));
            end
        end
        function display_iso(R)
            ISO = reshape(permute(R.S.ISO,[1 3 2]),[R.S.nt*R.S.nexp 3]);
            M1 = max(abs(ISO(:,1))); M2 = max(max(abs(ISO(:,2:3))));
            ISO(:,1) = ISO(:,1)*M2/M1/2;
            
            tidx = (1:R.S.nt*R.S.nexp)/R.S.nt + .5;
            cols = [1 1 1; .7 .7 1; 1 .7 .7];
            b = ones(1,R.S.nexp)*2; % code: 1=not selected, 2=rest, 3=stim
            b(R.S.stim)=3;
            b(~R.S.trials)=1;
            im = shiftdim(cols(b,:),-1);
            R.trimg = image([1 R.S.nexp],[-M2 M2]/3,im,'parent',R.grob.haiso); 
            axis(R.grob.haiso,'normal') % cancel the 'axis image' from function fn_imvalue
            line(tidx,ISO,'parent',R.grob.haiso,'hittest','on', ...
                'buttondownfcn',@(h,evnt)trial_toggle(R,R.grob.haiso));
        end
        function display_result(R,klist)
            if nargin<2, klist = 1:4; end
            
            % shortcuts
            doavgtrials = R.S.parameters.fastmain.resultavgtrials;
            
            % loop on displays
            for k=klist
                % display parameters: note that if only the current trial
                % is displayed, the flag about stim and/or rest has to be
                % ignored
                pars = R.options.resdisp(k);
                if strcmp(pars.trials,'current')
                    pars.trials = num2str(R.ktrial);
                    par.res = 'all'; 
                end
                
                % get results 
                datac = getresult(R.S,[R.kY R.kedg],pars);
                data = cat(1,datac{:})';
                if isempty(data)
                    ha = R.grob.hares(k); cla(ha)
                    if strcmp(pars.data,'volume')
                        title(ha,'no data available')
                    else
                        title(ha,'no result available')
                    end
                    continue
                end
                [n1 n2] = size(datac);
                % add volume results in a ROI if required
                if R.dispvolroi && strcmp(pars.data,'volume')
                    vol = fn_switch(pars.res, ...
                        'all',      mean(R.volumeroi,2), ...
                        'stim/rest',[R.volumeroi(:,1) R.volumeroi(:,2)], ...
                        'stim+rest',[R.volumeroi(:,1) R.volumeroi(:,2)], ...
                        'stim',     R.volumeroi(:,1), ...
                        'rest',     R.volumeroi(:,2));
                    datac = [datac vol]; %#ok<AGROW>
                end

                % display
                ha = R.grob.hares(k);
                % plot volume ROI results in gray if required
                if R.dispvolroi && strcmp(pars.data,'volume')
                    plot(1:nt,vol,'--','parent',ha,'color',[1 1 1]*.5)
                    hold(ha,'on')
                end
                % correction by edges interpolation spacing
                if R.dxcorr, data = data*R.dx; end
                % now plot vessel results
                hl=plot(data,'parent',ha,'linewidth',1.5);
                hold(ha,'off')
                axis(ha,'tight')
                % title
                title(ha,[pars.data ' - ' pars.res ' - ' pars.trials])
                % legend
                L = struct('x',datac, ...  % structure of same size as datac!
                    'data','','res','','trials','','name','');
                if strcmp(pars.data,'compare')
                    resname = R.S.parameters.resultnames;
                    for i=1:length(resname), [L(:,:,i).data] = deal(resname{i}); end
                end
                if strcmp(pars.trials,'odd+even')
                    [L(:,1,:).trials] = deal('odd');
                    [L(:,2,:).trials] = deal('even');
                end 
                if strcmp(pars.res,'stim+rest')
                    [L(1,:,:).res] = deal('stim');
                    [L(2,:,:).res] = deal('rest');
                end
                for i=1:numel(L)
                    tmp = {L(i).data L(i).res L(i).trials};
                    for j=1:3
                        a = tmp{j};
                        if isempty(a)
                        elseif isempty(L(i).name)
                            L(i).name = a;
                        else
                            L(i).name = [L(i).name '-' a];
                        end
                    end
                end
                if numel(L)>1, legend(hl,L.name,'location','southwest'), end 
                % callback
                if ~doavgtrials
                    set(hl,'hittest','on','buttondownfcn',@(h,evt)trial_toggle(R,ha))
                end
            end
        end
        function display_trialest(R)
            ha = R.grob.halin;
            ax = axis(ha(1));
            % volume
            data = getdataf(R.e,R.ktri);
            % flow estimation result
            res = {getJ(R.e,R.ktri),getV(R.e,R.ktri)};
            % shift pixels?
            npix = get(R.grob.hulinshift,'value')-1;
            if npix>10, npix=-(npix-10); end
            if npix
                np = R.e.np; nt = R.S.nt;
                np2 = np+npix;
                [ii jj] = ndgrid(1:np,1:nt);
                ii = fn_add(ii,-(0:nt-1)*npix); % this is the pixel shift
                ii = 1+mod(ii-1,np2);
                idx = ii + np2*(jj-1);
                tmp = [{data} res];
                for i=1:3
                    if ~isempty(tmp{i})
                        a = ones(np2,nt)*min(tmp{i}(:));
                        a(idx) = tmp{i};
                        tmp{i} = a;
                    end
                end
                data = tmp{1}; res = tmp(2:3);
            end
            % display
            if isempty(data)
                title(ha(1),'no data available'), cla(ha(1))
            else
                imagesc(data,'parent',ha(1))
            end
            if any(ax(1:2)~=[0 1]), fn_imvalue('chgx',ax(1:2),ha(1)), end
            if isempty(res{1})
                title(ha(2),'no result available'), cla(ha(2))
                title(ha(3),'no result available'), cla(ha(3))
            else
                imagesc(res{1},'parent',ha(2))
                imagesc(res{2},'parent',ha(3))
            end
            set(ha,'ydir','normal')
        end
    end
    
    % actions
    methods
        function trial_toggle(R,ha)
            p = get(ha,'currentpoint');
            itrial = round(p(1));
            if R.chgtrials && strcmp(get(gcf,'selectiontype'),'normal') ...
                    && ha==R.grob.haiso
                R.ktri = itrial; % no automatic display update!
                bval = ~R.S.trials(itrial);
                b0 = ones(1,R.S.nexp)*2; % code: 1=not selected, 2=rest, 3=stim
                b0(R.S.stim)=3;
                fn_buttonmotion(@trialtogglesub)
                saveproject(R.S) % save change in trials
                % update display
                display_result(R), drawnow
                display_trialest(R)
            else
                R.ktrial = itrial; % automatic display update
            end
            function trialtogglesub
                p = get(R.grob.haiso,'currentpoint');
                itrial2 = min(R.S.nexp,max(1,round(p(1))));
                idx = sort([itrial itrial2]);
                R.S.trials(idx(1):idx(2)) = bval;
                b = b0;
                b(~R.S.trials) = 1;
                cols = [1 1 1; .7 .7 1; 1 .7 .7];
                im = shiftdim(cols(b,:),-1);
                set(R.trimg,'cdata',im)
            end
        end
        function axesclick(R,ha)
            if ismember(ha,R.grob.halin(1:2))
                switch get(gcf,'SelectionType')
                    case 'open'
                        V = getV(R.e,R.ktri);
                        if ~any(V(:)), return, end
                        p = get(ha,'currentpoint'); p = p(1,1:2)';
                        [nx nt] = size(V);
                        t = round(p(1)); x = p(2);
                        if t<1 || t>nt || x<1 || x>nx, return, end
                        t0 = t; x0 = x; xa = []; xb = [];
                        while t>=1 && x>=1 && x<=nx
                            vxt = interp1(1:nx,V(:,t),x);
                            x = x + vxt;
                            t = t-1;
                            xa(end+1) = x; %#ok<AGROW>
                        end
                        t = t0; x = x0;
                        while t<=nt && x>=1 && x<=nx
                            vxt = interp1(1:nx,V(:,t),x);
                            x = x - vxt;
                            t = t+1;
                            xb(end+1) = x; %#ok<AGROW>
                        end
                        xdata = [xa(end:-1:1) x0 xb];
                        tdata = t0 + (-length(xa):length(xb));
                        hl(1) = line(tdata,xdata,'color','y','parent',R.grob.halin(1));
                        hl(2) = line(tdata,xdata,'color','y','parent',R.grob.halin(2));
                        set(hl,'hittest','on','buttondownfcn',@(h,evnt)trydelete(hl))
                end
            end
        end
    end
    
    % volume ROI
    methods
        function volumeROI(R)
            % load volume results
            fname = [R.S.savedir 'volume.mat'];
            if ~exist(fname,'file'), errordlg('no volume results'), return, end
            load(fname)
            hf = figure('numbertitle','off','name','define ROI...');
            [nx ny nt] = size(volstim); %#ok<NODEF>
            % draw ROI
            imagesc(volstim(:,:,1)')
            poly = fn_mouse('poly');
            % normalize volume results by average frame
            volall = (volstim+volrest)/2; %#ok<NODEF>
            avg = mean(volall,3);
            volstim = fn_mult(volstim,1./avg);
            volrest = fn_mult(volrest,1./avg);
            % get the signal
            mask = poly2mask(poly(2,:),poly(1,:),nx,ny);
            vol = cat(4,volstim,volrest);
            vol = reshape(vol,[nx*ny nt 2]);
            vol = squeeze(mean(vol(find(mask),:,:),1)); %#ok<FNDSB>
            % time interpolation to compensate for binning
            voltime = fn_bin(1:R.S.nt,volbin(2));
            vol = interp1(voltime,vol,1:R.S.nt,'linear','extrap');
            % display time courses
            plot(1:R.S.nt,vol), drawnow
            R.volumeroi = vol;
            % update display if necessary
            if R.dispvolroi, display_result(R), end
            % close figure
            close(hf)
        end
    end
    
    % misc
    methods
        function access(R) %#ok<MANU>
            keyboard
        end
    end
end

%---
function trydelete(hl)

delete(hl(ishandle(hl)))

end
    
    

