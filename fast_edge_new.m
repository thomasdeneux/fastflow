classdef fast_edge < hgsetget
   
    % PROPERTIES
    properties (SetAccess='private')
        % points
        snake
        halfd
        tubea
        tubeb
        dx = .8; % distance between points (for backward compatibility, default .8 value which was used until now)
        % parameters
        np
        nt
        nexp
    end
    properties
        % user
        active = true;
        flag = '';
        note = '';
    end
    properties (Dependent, SetAccess='private', Transient)
        % volume data
        data
    end
    properties
        % estimations
        flow = fast_flow.empty(1,0);
    end
    properties (Dependent, SetAccess='private', Transient)
        % parameters
        nest
    end
    properties
        % file info
        savedir
        filebase
        fakepar
        id
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
    
    % Constructor, copy
    methods
        function E = fast_edge(S,snake,halfd,dx,varargin)
            % function E = fast_edge()
            % function E = fast_edge(S,snake,halfd,dx[,'fake',fakepar])
            if nargin==0, return, end
            % checks
            if size(snake,1)~=2, error('snake should be a 2-rows array'), end
            if size(halfd,1)~=1, error('diameter should be a row vector'), end
            if size(snake,2)~=size(halfd,2)
                error('snake and diameter should have the same number of columns')
            end
            % set properties
            E.snake = snake;
            E.halfd = halfd;
            E.dx = dx;
            E.savedir = S.savedir;
            E.np   = size(snake,2);
            E.nt  = S.nt;
            E.nexp = S.nexp;
            E.id = fn_hash({snake,halfd},6);
            E.filebase = ['edge' E.id];
            if nargin>=5
                E.fakepar = struct('name',varargin{5},'value',varargin{6},'id',[]);
                E.fakepar.id = fn_hash(E.fakepar.value,6);
                E.filebase = [E.filebase '_' E.fakepar.name E.fakepar.id];
                fil = [E.fakepar.name E.fakepar.id '.xml'];
                if ~exist([E.savedir fil],'file')
                    disp(['write fake parameters in file ' fil])
                    fn_savexml([E.savedir fil],E.fakepar.value);
                end
            end
        end
        function updatesavedir(E,d)
            for i=1:length(E)
                E.savedir = d;
                updatesavedir(E.flow,d)
            end
        end
        function E = copy(E0) %#ok<MANU,STOUT>
            error('not implemented yet')
            if ~isscalar(E0) %#ok<UNRCH>
                if ndims(E0)>2, error('cannot copy ND arrays of edges'), end
                s = size(E0);
                E = fast_edge;
                E(s(1),s(2)) = fast_edge;
                for k=1:prod(s), E(k) = copy(E0(k)); end
                return
            end
            E = fast_edge;
            E.snake = E0.snake;
            E.halfd = E0.halfd;
            E.tubea = E0.tubea;
            E.active = E0.active;
            E.flag = E0.flag;
            E.savedir = E0.savedir;
            E.interppar = E0.interppar;
            E.filterpar = E0.filterpar;
            E.estpar = E0.estpar;
            E.estpar_local = E0.estpar_local;
            E.binpar = E0.binpar;
            E.savepar = E0.savepar;
        end
    end
    
    % Vessel: points and flags
    methods
        function[tubea tubeb] = gettube(E)
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
        function [col col0] = getcolor(E)
            col0 = fn_switch(E.flag,'','b', ...
                'good','y','medium',[1 .6 0],'bad','r');
            col = fn_switch(E.active,col0,'k');
        end
    end
    
    % Flow
    methods
        function n = get.nest(E)
            n = length(E.flow);
        end
        function addflow(E,filterpar,estpar)
            methods = regexp(estpar.method,'(radon|track)','tokens');
            for i=1:length(methods)
                method = methods{i}{1};
                par0 = feval(['fast_' method],'par'); % e.g. fast_track('par')
                estpar = fn_structmerge(par0,estpar,'skip');
                s = struct('METHOD',method, ...
                    'FILTER',{filterpar},'ESTPAR',estpar);
                idflow = fn_hash(s,3);
                iflow = find(strcmp({E.flow.id},idflow),1);
                if ~isempty(iflow)
                    % flow element already exists with these estimation parameters
                    continue
                end 
                fil = [method idflow '.xml'];
                if ~exist([E.savedir fil],'file')
                    disp(['write estimation parameters in file ' fil])
                    fn_savexml([E.savedir fil],s);
                end
                iflow = length(E.flow)+1;
                E.flow(iflow) = fast_flow(method,filterpar,estpar, ...
                    idflow,E.savedir,E.filebase);
            end
        end
    end
    
    % Data: file access
    methods
        function fname = mapfile(E,docreate)
            % function fname = mapfile(E,docreate)
            if nargin<2, docreate = false; end
            fil = [E.filebase '_data'];    
            fname = [E.savedir fil];
            if ~exist(fname,'file')
                if docreate
                    % create file if necessary
                    disp(['create file ' fil])
                    fid = fopen([E.savedir fil],'w');
                    fwrite(fid,zeros(fname,'single'),'single');
                    fclose(fid);
                else
                    fname = '';
                end
            end
        end
        function x = getdata(E,f,varargin)
            % function x = getdata(E,f[,whichtrials][,'avg|'sum'][,classname])
           
            % inputs
            ktrial = 0; sumflag = false; avgflag = false; classname = 'single';
            for i=1:length(varargin)
                a = varargin{i};
                if ischar(a)
                    switch a
                        case 'sum'
                            sumflag = true;
                        case 'avg'
                            sumflag = true;
                            avgflag = true;
                        otherwise
                            classname = a;
                    end
                else
                    ktrial = a;
                end
            end
            readformat = ['single=>' classname];
            
            % file name and data size
            fname = mapfile(E,false);
            if isempty(fname), x=[]; return, end
            siz = [E.np E.nt E.nexp];
            
            % read file: in a block or trial by trial
            fid = fopen(fname,'r');
            if all(ktrial==0) % not really faster to read a single block
                x = reshape(fread(fid,prod(siz),readformat),siz);
                if sumflag, x = sum(x,3); end
                if avgflag, x = x / siz(3); end
            elseif ~sumflag
                x = zeros(siz(1),siz(2),length(ktrial),classname);
                for i=1:length(ktrial)
                    fseek(fid,(ktrial(i)-1)*prod(siz(1:2))*4,'bof');
                    x(:,:,i) = fread(fid,siz(1:2),readformat);
                end
            else
                x = zeros(siz(1:2),classname);
                for k=ktrial
                    fseek(fid,(k-1)*prod(siz(1:2))*4,'bof');
                    x = x + fread(fid,siz(1:2),readformat);
                end
                if avgflag, x = x / length(ktrial); end
            end
            fclose(fid);
        end
        function setdata(E,ktrial,x)
            % function setdata(E,ktrial,x)
         
            % file name and size
            fname = mapfile(E,true);
            siz = [E.np E.nt E.nexp];
            
            % write in file
            fid = fopen(fname,'r+');
            fseek(fid,(ktrial-1)*prod(siz(1:2))*4,'bof');
            fwrite(fid,x,'single');
            fclose(fid);
        end
        function b = hasdata(E,ktrial)
            % function b = hasdata(E,ktrial)
            
            % multiple E?
            if length(E)>1
                b = false(length(E),length(ktrial));
                for i=1:length(E)
                    b(i,:) = hasdata(E(i),ktrial);
                end
                if isvector(b), b = b(:)'; end
                return
            end
            
            % singleton E
            b = false(1,length(ktrial));
            % file name
            fname = mapfile(E);
            if isempty(fname), return, end
            % read file
            fid = fopen(fname,'r+');
            for i=1:length(ktrial)
                fseek(fid,(ktrial-1)*prod(siz(1:2))*4,'bof');
                b(i) = logical(fread(fid,1,'single'));
            end
            fclose(fid);
        end
        function x = get.data(E)
            x = getdata(E);
        end
    end
    
    % Data: filtering
    methods
        function y = getdataf(E,k)
            y = getdata(E,k,'double');
            for i=1:length(E.filterpar)
                eval(E.filterpar{i});
            end
        end
    end
    
    % Backward compatibility
    methods (Static)
        function E = loadobj(x)
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
                    'savedir',[x.savedir 'new2/']);
                isfake = isfield(x.interppar,'fakepar');
                if isfake, fakearg={'fake',x.interppar.fakepar}; else fakearg={}; end
                % set edge
                E = fast_edge(s,x.snake,x.halfd,.8,fakearg{:});
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
                filedatanew = [E.savedir E.filebase '_data'];
                if exist(filedataold,'file')
                    disp([' rename ' filedataold ' to ' filedatanew])
                    if exist(E.filebase,'file')
                        fprintf('\b (error: file exists)\n')
                    else
                        copyfile(filedataold,filedatanew)
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
                            'FILTER',{sxml.FILTER},'ESTPAR',estpark);
                        idflownew = fn_hash(s,3);
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
                    %delete(fnameold)
                end
                disp('finished, now please remove old xml files manually')
                cd(swd)
                disp(' ')
            else
                E = x;
            end
        end
    end
    
    
end



