classdef fast_edge < fast_data
    % function E = fast_edge()
    % function E = fast_edge(S,snake,halfd,dx)
    % function E = fast_edge(S,'fake',fakepar)
    
    % PROPERTIES
    properties (SetAccess='private')
        % points
        snake
        halfd
        tubea
        tubeb
        dx = .8; % distance between points (for backward compatibility, default .8 value which was used until now)
    end
    properties
        % estimations
        flow = fast_flow.empty(1,0);
        % section estimations
        section = fast_section.empty(1,0);
    end
    properties (Dependent, SetAccess='private', Transient)
        % data
        data
        % parameters
        nest
    end
    properties (SetAccess='private')
        % additional file info
        interpflag = '';
        filepar
        fakepar
    end
    properties (Access='private')
        % for backward compatibility only
        version = 2;
        interppar
        filterpar
        estpar
        binpar
        savepar
    end
    properties (SetAccess='private')
        flagnumber = 1;
    end
    properties (Dependent)
        flag
    end
    
    % Constructor, copy
    methods
        function E = fast_edge(S,varargin)
            % function E = fast_edge()
            % function E = fast_edge(S,snake,halfd,dx)
            % function E = fast_edge(S,'fake',fakepar)
            % function E = fast_edge(S,'scan',fname)
            if nargin==0, arg={}; else arg={S}; end
            E = E@fast_data(arg{:});
            if nargin==0, return, end
            % set properties
            if isnumeric(varargin{1})
                E.interpflag = 'movie';
                [E.snake E.halfd E.dx] = deal(varargin{:});
                % checks
                if size(E.snake,1)~=2, error('snake should be a 2-rows array'), end
                if size(E.halfd,1)~=1, error('diameter should be a row vector'), end
                if size(E.snake,2)~=size(E.halfd,2)
                    error('snake and diameter should have the same number of columns')
                end
                E.np = size(E.snake,2);            
                E.nt = S.nt;
                E.nexp = S.nexp;
                E.dx = 1;
                E.id = fn_hash({E.snake,E.halfd},6);
                E.filebase = ['edge' E.id];
            elseif strcmp(varargin{1},'scan')
                E.interpflag = 'scan';
                files = varargin{2};
                oksame = ~any(diff(double(files),1,1),1);
                filebase = files(1,:);
                filebase(~oksame) = 'x';
                E.filepar = struct('filebase',filebase,'files',files);
                a = fn_readimg([S.datadir deblank(files(1,:))]);
                [E.np E.nt dum] = size(a);
                E.nexp = size(files,1);
                E.dx = 1;
                E.id = fn_hash(E.fakepar,6);
                E.filebase = ['scan' E.id];
                fil = [E.filebase '.xml'];
                if ~exist([E.savedir fil],'file')
                    disp(['write fake parameters in file ' fil])
                    fn_savexml([E.savedir fil],E.fakepar);
                end
            elseif strcmp(varargin{1},'fake')
                E.interpflag = 'fake';
                E.fakepar = varargin{2};
                E.np = E.fakepar.nx;
                E.nt = E.fakepar.nt;
                E.nexp = E.fakepar.nexp;
                E.id = fn_hash(E.fakepar,6);
                E.filebase = ['fake' E.id];
                fil = [E.filebase '.xml'];
                if ~exist([E.savedir fil],'file')
                    disp(['write fake parameters in file ' fil])
                    fn_savexml([E.savedir fil],E.fakepar);
                end
            else
                error argument
            end
            E.algo = 'volume';
        end
        function E = copy(E0) %#ok<MANU,STOUT>
            error('not implemented yet')
        end
    end
    
    % Vessel: points, interpolation array and flags
    methods
        function[tubea tubeb] = gettube(E)
            if isempty(E.snake), tubea=[]; tubeb=[]; return, end
            dsnake = gradient(E.snake);
            dnorm  = sqrt(sum(dsnake.^2));
            snaken = [dsnake(2,:)./dnorm; -dsnake(1,:)./dnorm];
            hdiam = repmat(E.halfd,2,1);
            tubea = E.snake + snaken.*hdiam;
            tubeb = E.snake - snaken.*hdiam;
        end
        function x = get.tubea(E)
            if isempty(E.tubea)
                [E.tubea E.tubeb] = gettube(E);
            end
            x = E.tubea;
        end
        function x = get.tubeb(E)
            if isempty(E.tubeb)
                [E.tubea E.tubeb] = gettube(E);
            end
            x = E.tubeb;
        end
    end
    
    % Vessels: interpolation
    methods
        function y = interptube(E,Y)
            % function y = interptube(E,Y)
            %---
            % interpolate data from a nx*ny*nt array Y using the vessel
            % tube, resulting in an np*nt array
            y = fast_interptube(Y,E);
        end
        function a = sparsemerging(E,nx,ny)
            % function a = sparsemerging(E,nx,ny)
            %---
            % returns a sparse array to merge vessel data into a nx*ny
            % image
            ne = length(E);
            a = cell(1,ne);
            for i=1:ne
                e = E(i);
                tube = double([e.tubea e.tubeb(:,end:-1:1)]);
                mask = poly2mask(tube(2,:),tube(1,:),nx,ny);
                [ii jj] = find(mask); % column vectors
                dd = fn_add(ii,-e.snake(1,:)).^2 + fn_add(jj,-e.snake(2,:)).^2;
                [dum kk] = min(dd,[],2);
                a{i} = sparse(ii+nx*(jj-1),kk,1,nx*ny,e.np);
            end
            a = [a{:}];
            a1 = sum(a,2)'; % several vessels can contribute to same pixel
            for i=find(a1>1), a(i,:) = a(i,:)/a1(i); end
        end
    end
    
    % Flag & color
    methods (Static)
        function s = flaglist
            cols = [0 0 1; 1 1 0; 1 .5 0; 1 0 0; ...
                .5 0 0; .5 .5 0; 0 1 0; 0 .5 0; ...
                0 1 1; 0 .5 .5; 0 0 .5; ...
                1 0 1; .5 0 .5];
            ncol = size(cols,1);
            flags = {'' 'good' 'medium' 'bad'};
            nflag = length(flags);
            flags(nflag+1:ncol) = num2cell(char('a'+(nflag+1:ncol)-1));
            s = struct('flag',flags,'color',num2cell(cols,2)');
        end
    end
    methods
        function k = get.flagnumber(E)
            k = E.flagnumber;
            n = length(E.S.flaglist.alllists); % number of flag lists
            k(end+1:n) = 1;
        end
        function flag = getflag(E)
            if ~isa(E.S,'fastflow'), flag=[]; return, end
            klist = E.S.flaglist.k; % index of current flag list
            flaglist = getflaglist(E.S); % current flag list
            flag = flaglist(E.flagnumber(klist));
        end
        function k = getflagnumber(E)
            klist = E(1).S.flaglist.k; % index of current flag list
            k = cat(1,E.flagnumber);
            k = k(:,klist)';
        end
        function setflagnumber(E,k)
            % function setflagnumber(E,k)
            %---
            % modifies the flag for the current flag list (index given by
            % E.S.flaglist.k)
            klist = E.S.flaglist.k; % index of current flag list
            nlist = length(E.flagnumber);
            E.flagnumber(klist) = k;
            E.flagnumber(nlist+1:klist)=1;
        end
        function flag = get.flag(E)
            flag = getflag(E);
            if isempty(flag)
                flag = '';
            else
                flag = flag.flag;
            end
        end 
        function [col col0] = getcolor(E)
            flag = getflag(E);
            if isempty(flag)
                col0 = [0 0 0];
            else
                col0 = flag.color;
            end
            col = fn_switch(E.active,col0,[0 0 0]);
        end
    end
    
    % Flow
    methods
        function n = get.nest(E)
            n = length(E.flow);
        end
        function iflows = addflow(E,varargin)
            % function iflows = addflow(E,method,filterpar,estpar)
            % function iflows = addflow(E,filterpar,estpar)
            % function iflows = addflow(E,idflow)
            %---
            % in the second case, estpar must have a field 'method'
            % in the third case, no new flow is created
            switch nargin
                case 2
                    idflow = varargin{1};
                    iflows = find(strcmp({E.flow.id},idflow),1);
                    return
                case 3
                    [filterpar estpar] = deal(varargin{:});
                    methods = regexp(estpar.method,'(radon|track)','tokens');
                case 4
                    [method filterpar estpar] = deal(varargin{:});
                    methods = {{method}};
                otherwise
                    error('too many arguments')
            end
            iflows = zeros(1,length(methods));
            for i=1:length(methods)
                method = methods{i}{1};
                par0 = feval(['fast_' method],'par'); % e.g. fast_track('par')
                par1 = fn_structmerge(par0,estpar,'skip'); 
                [idflow s] = fast_flow.getflowid(method,filterpar,par1); 
                iflow = find(strcmp({E.flow.id},idflow),1);
                if ~isempty(iflow)
                    % flow element already exists with these estimation parameters
                    iflows(i) = iflow;
                    continue
                end 
                fil = [method idflow '.xml'];
                 if ~exist([E.savedir fil],'file')
                    disp(['write estimation parameters in file ' fil])
                    fn_savexml([E.savedir fil],s);
                 end
                iflow = length(E.flow)+1;                
                E.flow(iflow) = fast_flow(E,method,filterpar,par1,idflow); 
                iflows(i) = iflow;
            end
            if nargout==0, clear iflows, end
        end
    end
    
    % Vessel section
    methods
        function isection = addsection(E,a)
            % function isection = addsection(E,sectionpar)
            % function isection = addsection(E,idsection)
            %---
            % in the second case, no new fast_section object is created
            if ischar(a)
                idsection = a;
                isection = find(strcmp({E.section.id},idsection),1);
                return
            end
            sectionpar = a;
            par0 = fast_vesselsection('par');
            par1 = fn_structmerge(par0,sectionpar,'skip');
            idsection = fast_section.getsectionid(par1); 
            isection = find(strcmp({E.section.id},idsection),1);
            if ~isempty(isection)
                return
            end
            fil = ['section' idsection '.xml'];
            if ~exist([E.savedir fil],'file')
                disp(['write section parameters in file ' fil])
                fn_savexml([E.savedir fil],par1);
            end
            isection = length(E.section)+1;
            E.section(isection) = fast_section(E,par1,idsection);
            if nargout==0, clear isection, end
        end
    end
    
    % Data: file access
    methods
        function x = getdata(E,a)
            % function x = getdata(E,ktrials)
            % function x = getdata(E,resdisp)
            if isstruct(a)
                resdisp = a;
            else
                resdisp = struct('data','volume','trial',a);
            end
            x = mapgetresult(E,resdisp);
        end
        function pos = getrbcpos(E,ktrial)
            % function x = getrbcpos(E,ktrial)
            if ~isnumeric(ktrial) || ~isscalar(ktrial)
                error('getrbcpos(E,ktrial) used with non-numeric or non-scalar ktrial')
            end
            ok = false; k=0; pos = {};
            while ~ok
                k = k+1;
                block = mapread(E,ktrial,sprintf('_rbcpos%i',k));
                pos{k} = block; %#ok<AGROW>
                ok = ~block(E.np,1);
            end
            klast = find(block(:,1),1,'last');
            pos{k} = block(1:klast,:);
            pos = cat(1,pos{:});
        end
        function setdata(E,ktrial,x,pos)
            % function setdata(E,ktrial,x[,pos])
            %---
            % x is an np x nt array
            % pos is an nRBC x nt array, used with fake data (it must
            % comply nRBC<=np)
            mapwrite(E,ktrial,x)
            if isfield(E.tmpdata,'current_volume') && E.tmpdata.current_volume.ktrial==ktrial
                E.tmpdata = rmfield(E.tmpdata,'current_volume');
            end 
            if nargin>=4
                [nRBC nt1] = size(pos);
                if nt1~=E.nt
                    error('size of RBC positions array must be nRBC x nt')
                end
                nblocks = ceil((nRBC+1)/E.np);
                pos(nblocks*E.np,1)=0;
                for k=1:nblocks
                    block = pos((k-1)*E.np+(1:E.np),:);
                    mapwrite(E,ktrial,block,sprintf('_rbcpos%i',k))
                end
            end
        end
        function setdatatmp(E,ktrial,x)
            % function setdatatmp(E,ktrial,x)
            settmpdata(E,ktrial,'volume',x)
        end
        function b = hasdata(E,ktrial,doforce)
            % function b = hasdata(E,ktrial[,doforce])
            if nargin<3, doforce=false; end
            b = maphascontent(E,ktrial,doforce);
        end
        function x = get.data(E)
            x = mapread(E,[]);
        end
        function x = getresult(E,resdisp)
            % function x = getresult(E,resdisp)
            %---
            % finds the appropriate object (fast_edge or fast_data) to get
            % result from before calling fast_data.mapgetresult on this object
            %
            % resdisp.algo can be:
            % - '' (when resdisp.data is 'volume')
            % - 'result' for the first result
            % - method (e.g. 'track') for the first result with this method
            % - methodcode (e.g. 'trackABCD') for the restul with this algo
            % - add 'lines' (e.g. 'track_lines') to get the lines result
            %
            % see also fast_flow.getresult
            
            % multiple edges -> return cell array
            if ~isscalar(E)
                x = cell(size(E));
                for i=1:numel(E)
                    x{i} = getresult(E(i),resdisp);
                end
                return
            end
        
            
            % volume
            switch resdisp.data
                case {'volume' 'flow'}
                    % get volume data from this fast_edge object
                    x = mapgetresult(E,resdisp);
                case {'volumef','result','lines'}
                    % get filtered volume data or flow result from the
                    % appropriate fast_flow child 
                    if isempty(E.flow), x=[]; return, end
                    % which algo
                    if isfield(resdisp,'algo'), algo=resdisp.algo; else algo=1; end
                    % which fast_flow child to read result from?
                    if isnumeric(algo)
                        kflow = algo;
                    else
                        % decompose in method and id
                        tmp = regexp(algo,'^([a-z]*)([A-P]*)$','tokens');
                        if isempty(tmp), error('wrong algo specification'), end
                        [method idflow] = deal(tmp{1}{:});
                        if isempty(idflow)
                            methods = {E.flow.method};
                            kflow = find(strcmp(methods,method),1);
                            if isempty(kflow), x=[]; return, end % no result with this method
                        else
                            algos = {E.flow.algo};
                            kflow = find(strcmp(algos,algo),1);
                            if isempty(kflow), error('no result ''%s''',algo), end
                        end
                    end
                    % get result!
                    x = mapgetresult(E.flow(kflow),resdisp);
                case {'section','sigma'}
                    % get section data or estimated width from the
                    % appropriate fast_section child
                    if isempty(E.section), x=[]; return, end
                    x = mapgetresult(E.section(1),resdisp);
                otherwise
                    error('unknown data flag ''%s''',resdisp.data)
            end
        end
        function chgnexp(E,n)
            % function chgnexp(E,n)
            
            % multiple E
            if ~isscalar(E), for i=1:length(E), chgnexp(E(i),n), end, return, end
            if isempty(E.fakepar)
                error('property ''nexp'' can be changed only in fake mode')
            end
            % check
            if n<E.nexp
                b = hasresult(E.flow,n+1:E.nexp,true);
                nall = n+find(all(b,1),1,'last');
                nany = n+find(any(b,1),1,'last');
                if isempty(nany)
                    % no result at all: ok to remove trials
                    answ = 'Continue';
                elseif isempty(nall)
                    if nany==E.nexp
                        answ = questdlg({'Some results exist until the last trial,' ...
                            'do you want to loose them?'},'confirmation','Continue (yes)','Cancel','Cancel');
                    else
                        answ = questdlg({['Some results exist until trial ' num2str(nany) ','] ...
                            ['do you prefer to reduce at ' num2str(nany) ' instead of ' num2str(n) '?']}, ...
                            'confirmation',num2str(n),num2str(nany),'Cancel','Yes',num2str(nany));
                    end
                else
                    if nall==E.nexp
                        answ = questdlg({'All results exist until the end for some algorithms,' ...
                            'do you want to loose them?'},'confirmation','Continue (yes)','Cancel','Cancel');
                    elseif nany>nall && nany==E.nexp
                        answ = questdlg({['All results exist until trial ' num2str(nall) ', '], ...
                            'some results exist until the last trial,', ...
                            ['do you prefer to reduce at ' num2str(nall) ' instead of ' num2str(n) ', or cancel?']}, ...
                            'confirmation',num2str(n),num2str(nall),'Cancel','Cancel');
                    elseif nany>nall
                        twonum = [num2str(nall) ' or ' num2str(nany)];
                        answ = questdlg({['All results exist until trial ' num2str(nall) ', '], ...
                            ['some results exist until trial ' num2str(nany) ','], ...
                            ['do you prefer to reduce at ' num2str(nall) ' or ' num2str(nany) ' instead of ' num2str(n) '?']}, ...
                            'confirmation',num2str(n),twonum,'Cancel',twonum);
                        if strcmp(answ,twonum)
                            answ = questdlg('Please choose:', ...
                                'confirmation',num2str(nany),num2str(nall),'Cancel',num2str(nany));
                        end
                    else
                        answ = questdlg({['All results exist until trial ' num2str(nall) ', '], ...
                            ['do you prefer to reduce at ' num2str(nall) ' instead of ' num2str(n) '?']}, ...
                            'confirmation',num2str(n),num2str(nall),'Cancel',num2str(nall));
                    end
                end
                if strcmp(answ,'Cancel')
                    return
                elseif findstr(answ,'Continue')
                    % confirmation: n is unchanged
                else
                    n = str2double(answ);
                    if isnan(n), error programming, end
                end
            end
            mapchgnexp(E,n)
            for i=1:length(E.flow), mapchgnexp(E.flow(i),n), end
        end
    end
    
    % Backward compatibility
    methods (Static)
        function E = loadobj(E)
            if isempty(E.interpflag)
                if isempty(E.S)
                    E.interpflag = '';
                elseif ~isempty(E.fakepar)
                    E.interpflag = 'fake';
                else
                    E.interpflag = 'movie';
                end
            end
        end
        function E = loadoldversion(x)
            % function E = loadoldversion(x)
            if isstruct(x)
                % load from previous version
                files = x.interppar.files;
                if isscalar(files)
                    n = files;
                elseif isvector(files) && islogical(files)
                    n = sum(files);
                else
                    n = size(files,1);
                end
                s = struct('nt',x.interppar.nfr,'nexp', n, ...
                    'savedir',x.savedir);
                isfake = isfield(x.interppar,'fakepar');
                if isfake, fakearg={'fake',x.interppar.fakepar}; else fakearg={}; end
                % set edge
                E = fast_edge(s,x.snake,x.halfd,.8,fakearg{:});
                if isfield(x,'flag'), E.flag   = x.flag; end
                if isfield(x,'active'), E.active = x.active; end
                % set flow subelements
                addflow(E,x.filterpar,x.estpar);
                % file transfer
                disp('updating files for new fast_edge object:')
                swd = pwd;
                cd(x.savedir)
                % (easy for edge data)
                nhash = fn_switch(isfake,6,3);
                id_data = fn_hash(x.interppar,nhash);
                filedataold = [E.filebase '_data' id_data];
                filedatanew = [E.savedir E.filebase '_volume'];
                if exist(filedataold,'file')
                    disp([' rename ' filedataold ' to ' filedatanew])
                    if exist(E.filebase,'file')
                        fprintf('\b (error: file exists)\n')
                    else
                        movefile(filedataold,filedatanew)
                    end
                else
                    disp(' no data file to rename')
                end
                % (difficult for estimation results data!)
                d = dir([E.filebase '_flow*']);
                filelist = {d.name};
                for k=1:length(filelist)
                    % interpret file name according to pattern
                    fnameold = filelist{k};
                    tmp = regexp(fnameold,['^' E.filebase '_flow' id_data '(.*)$'],'tokens');
                    if isempty(tmp)
                        disp([' error with file ' fnameold ': name does not match pattern'])
                        continue
                    end
                    id_flow = tmp{1}{1};
                    % read xml file
                    filexml = ['flow' id_flow '.xml'];
                    if ~exist(filexml,'file')
                        disp([' error with file ' fnameold ': cannot find xml file'])
                        continue
                    end
                    sxml = fn_readxml(filexml); % fields FILTER, EST, SAVE, BIN
                    if ~isequal(sxml.BIN,struct('xbin',1,'tbin',1))
                        disp([' error with file ' fnameold ': cannot deal with binning'])
                        continue
                    end
                    % add flow elements for this estimation result
                    addflow(E,sxml.FILTER,sxml.EST)
                    % copy old file into one or several new files
                    disp([' copy from file ' fnameold])
                    fin = fopen(fnameold,'r');
                    siz = [E.np E.nt E.nexp];
                    for i=1:length(sxml.SAVE)
                        % read
                        f = sxml.SAVE{i};
                        dat = reshape(fread(fin,prod(siz),'single=>single'),siz);
                        % new file name
                        method = fn_switch(f,'V','track','J','track', ...
                            'G','gabor','R','radon');
                        par0 = feval(['fast_' method],'par'); % e.g. fast_track('par')
                        estpark  = fn_structmerge(par0,sxml.EST,'skip');
                        s = struct('METHOD',method, ...
                            'FILTER',{sxml.FILTER(:)},'ESTPAR',estpark);
                        idflownew = fn_hash(s,4);
                        fnamenew = [E.savedir E.filebase '_' method idflownew];                      
                        if f=='J', fnamenew = [fnamenew '_lines']; end %#ok<AGROW>
                        % copy
                        disp([' -> in file ' fnamenew])
                        fout = fopen(fnamenew,'w');
                        fwrite(fout,dat,'single');
                        fclose(fout);
                    end
                    fclose(fin);
                    disp([' delete file ' fnameold])
                    delete(fnameold)
                end
                disp('finished, now please remove old xml files manually')
                cd(swd)
                disp(' ')
            else
                E = x;
            end
        end
    end
    
    % Misc
    methods
        function chgfield(E,fieldname,value)
            % function chgfield(E,fieldname,value)
            %---
            % applies both for fast_edge object(s) E and for fast_flow
            % children
            if isempty(E), return, end
            try
                chgfield@fast_data(E,fieldname,value)
                chgfield([E.flow],fieldname,value)
                chgfield([E.section],fieldname,value)
            catch
                [E.(fieldname)] = deal(value);
            end
        end
        function hidetmpdata(E,flag)
            % function hidetmpdata(E,'hide|show|remove')
            %---
            % applies both for fast_edge object(s) E and for fast_flow
            % children
           if isempty(E), return, end
            hidetmpdata@fast_data(E,flag)
            hidetmpdata([E.flow],flag)
            hidetmpdata([E.section],flag)
        end
        function recover(e,docleanfirst)
            % function recover(e[,docleanfirst])
            %---
            % updates the 'fast_flow' elements of e by checking the files
            % in the saving directory, put the right identification keys,
            % ...
            % if docleanfirst is set and is equal to true, deletes all the
            % existing 'fast_flow' elements first, and rebuild them all
            % based on existing files (this means that estimation
            % parameters for estimations which have not started yet will be
            % lost)
            if nargin<2, docleanfirst=false; end
            if docleanfirst
                [e.flow] = deal(fast_flow.empty);
            end
            swd = pwd;
            cd(e(1).savedir)
            for i=1:length(e)
                E = e(i);
                d = dir([E.filebase '_*']);
                filelist = {d.name};
                for k=1:length(filelist)
                    fname = filelist{k};
                    % interpret file name according to pattern
                    tmp = regexp(fname,['^' E.filebase '_([a-z])*([A-P])*$'],'tokens');
                    if isempty(tmp) || strcmp(tmp{1}{1},'volume')
                        continue
                    end
                    [method id_flow] = deal(tmp{1}{:});
                    % read xml file
                    filexml = [method id_flow '.xml'];
                    if ~exist(filexml,'file')
                        % create a new xml file
                        disp([' error with file ' fname ': cannot find xml file'])
                        continue
                    end
                    sxml = fn_readxml(filexml); % fields FILTER, EST, SAVE, BIN
                    % add flow elements for this estimation result
                    sxml.ESTPAR.method = method;
                    iflow = addflow(E,sxml.FILTER,sxml.ESTPAR);
                    % rename file if necessary
                    id_flow2 = E.flow(iflow).id;
                    if ~strcmp(id_flow,id_flow2)
                        % ID has been changed; rename file, and xml file
                        % should be trashed (but not now yet)
                        fname2 = strrep(fname,id_flow,id_flow2);
                        disp(['rename ' fname ' to ' fname2])
                        movefile(fname,fname2)
                        if exist([fname '_lines'],'file')
                            disp(['rename ' fname '_lines to ' fname2 '_lines'])
                            movefile([fname '_lines'],[fname2 '_lines'])
                        end
                    else
                        % mark the xml file as recently modified
                        system(['touch ' filexml]);
                    end
                end
            end
            cd(swd)
            disp('update finished, please remove old xml files manually')
        end
        function access(E) %#ok<MANU>
            keyboard
        end
        function eval(E,command) %#ok<MANU>
            eval(command)
        end
    end
    
end



