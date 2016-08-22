classdef fast_edge < hgsetget
   
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
        % user
        active = true;
        flag = '';
        note = '';
        % estimation parameters
        savedir
        interppar = struct('files',[],'nfr',0);
        filterpar = {};
        estpar
        estpar_local    % values in estpar_local are not used in the estimation but remind about parameters specific to that vessel (such as speed initialization)
        binpar  % binning of estimation results
        savepar % sizes of the estimation results (function of sizes and binpar)
    end
    properties (Dependent, SetAccess='private', Transient)
        np
        nfr
        nexp
    end
    properties (Dependent, SetAccess='private', Transient)
        data
        V
        J
        G
        R
    end
    properties (SetAccess='private', Transient)
        id_edge
        id_data
        id_algo
    end        
    
    % Constructor, copy
    methods
        function E = fast_edge(varargin)
            % function E = fast_edge(varargin)
            % function E = fast_edge(snake,halfd,dx,varargin)
            if nargin>0 && isnumeric(varargin{1})
                setsnake(E,varargin{1:3})
                varargin = varargin(4:end);
            end
            if ~isempty(varargin), set(E,varargin{:}), end
        end
        function E = copy(E0)
            if ~isscalar(E0)
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
    
    % GET/SET
    methods
        function setsnake(E,x,y,dx)
            if size(x,1)~=2, error('snake should be a 2-rows array'), end
            if size(y,1)~=1, error('diameter should be a row vector'), end
            if size(x,2)~=size(y,2)
                error('snake and diameter should have the same number of columns')
            end
            E.snake = x;
            E.halfd = y;
            E.dx = dx;
            E.id_edge = [];
        end
        function set.interppar(E,s)
            E.interppar = s;
            E.id_data = [];
        end
        function set.filterpar(E,str)
            if isempty(str)
                E.filterpar = {};
            else
                E.filterpar = cellstr(str);
            end
            E.id_algo = [];
        end
        function set.estpar(E,s)
            E.estpar = s;
            E.id_algo = [];
        end
        function set.binpar(E,s)
            E.binpar = s;
            E.id_algo = [];
        end
        function set.savepar(E,s)
            E.savepar = s;
            E.id_algo = [];
        end
        function n = get.np(E)
            n = size(E.snake,2);
        end
        function n = get.nfr(E)
            n = E.interppar.nfr;
        end
        function n = get.nexp(E)
            % possibility that interppar.files be either the name of files,
            % or a vector of booleans indicating which files are including,
            % or even just the number of files to use 
            files = E.interppar.files;
            if isscalar(files)
                n = files;
            elseif isvector(files) && islogical(files)
                n = sum(files);
            else
                n = size(files,1);
            end
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
        function h = get.id_edge(E)
            if isempty(E.id_edge)
               E.id_edge = fn_hash({E.snake,E.halfd},6);
            end
            h = E.id_edge;
        end
        function h = get.id_data(E)
            if isempty(E.id_data)
                nhash = fn_switch(isfield(E.interppar,'fakepar'),6,3);
                E.id_data = fn_hash(E.interppar,nhash);
                fil = ['data' E.id_data '.xml'];
                if exist(E.savedir,'dir') && ~exist([E.savedir fil],'file')
                    disp(['write interpolation parameters in file ' fil])
                    fn_savexml([E.savedir fil],E.interppar);
                end
            end
            h = E.id_data;
        end
        function h = get.id_algo(E)
            if isempty(E.id_algo)
                if isempty(E.savepar), names=[]; else names=fieldnames(E.savepar); end
                s = struct('FILTER',{E.filterpar},'EST',{E.estpar}, ...
                    'SAVE',{names},'BIN',{E.binpar});
                E.id_algo = fn_hash(s,3);
                fil = ['flow' E.id_algo '.xml'];
                if exist(E.savedir,'dir') && ~exist([E.savedir fil],'file')
                    disp(['write estimation parameters in file ' fil])
                    fn_savexml([E.savedir fil],s);
                end
            end
            h = E.id_algo;
        end
    end
    
    % file access
    methods
        function [fname siz offset] = mapfile(E,f,varargin)
            createflag = fn_flags('create',varargin);
            % function m = memmap(E,fieldname)
            siz = [E.np E.nfr E.nexp];
            % checks
            if any(siz==0) || isempty(E.interppar) || ~checksavedir(E) ...
                    || (~strcmp(f,'data') && (isempty(E.filterpar) ...
                    || isempty(E.estpar) || isempty(E.binpar) || isempty(E.savepar)))
                fname = []; siz = []; offset = [];
                return
            end
            if ~strcmp(f,'data') && ~isfield(E.savepar,f)
                fprintf('no estimation result named ''%s''\n',f)
                fname = []; siz = []; offset = [];
                return
            end
            % file name
            if strcmp(f,'data')
                fil = ['edge' E.id_edge '_data' E.id_data];
            else
                fil = ['edge' E.id_edge '_flow' E.id_data E.id_algo];
            end
            fname = fullfile(E.savedir,fil);
            % data size
            if strcmp(f,'data'), s = struct('data',siz); else s = E.savepar; end
            siz = s.(f);
            F = fieldnames(s);
            % create file if necessary
            if ~exist(fname,'file')
                if ~createflag, fname = []; siz = []; offset = []; return, end
                disp(['create file ' fil])
                fid = fopen(fname,'w');
                for i=1:length(F)
                    fwrite(fid,zeros(s.(F{i}),'single'),'single');
                end
                fclose(fid);
            end
            % offset
            offset = 0;
            for i=1:length(F)
                f1 = F{i};
                if strcmp(f1,f), break, end
                offset = offset + prod(s.(f1))*4;
            end
        end
        function x = mapread(E,f,varargin)
            % function x = mapread(E,f[,whichtrials][,'avg|'sum'][,classname])
           
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
            % data size
            [fname siz offset] = mapfile(E,f);
            if isempty(fname), x=[]; return, end
            % read file: in a block or trial by trial
            fid = fopen(fname,'r');
            if all(ktrial==0) % not really faster to read a single block
                fseek(fid,offset,'bof');
                x = reshape(fread(fid,prod(siz),readformat),siz);
                if sumflag, x = sum(x,3); end
                if avgflag, x = x / siz(3); end
            elseif ~sumflag
                x = zeros(siz(1),siz(2),length(ktrial),classname);
                for i=1:length(ktrial)
                    fseek(fid,offset + (ktrial(i)-1)*prod(siz(1:2))*4,'bof');
                    x(:,:,i) = fread(fid,siz(1:2),readformat);
                end
            else
                x = zeros(siz(1:2),classname);
                for k=ktrial
                    fseek(fid,offset + (k-1)*prod(siz(1:2))*4,'bof');
                    x = x + fread(fid,siz(1:2),readformat);
                end
                if avgflag, x = x / length(ktrial); end
            end
            fclose(fid);
        end
        function mapwrite(E,f,ktrial,x)
            % function mapwrite(E,f,ktrial,x)
         
            % checks
            [fname siz offset] = mapfile(E,f,'create');
            if isempty(fname)
                error('cannot write data: some properties might not be set correctly')
            end
            if size(x)~=siz(1:2), error('dimension mismatch'), end
            % write in file
            fid = fopen(fname,'r+');
            fseek(fid,offset + (ktrial-1)*prod(siz(1:2))*4,'bof');
            fwrite(fid,x,'single');
            fclose(fid);
        end
        function b = mapnonzero(E,f,ktrial)
            % function b = mapnonzero(E,f,ktrial)
            
            % multiple E?
            if length(E)>1
                b = false(length(E),length(ktrial));
                for i=1:length(E)
                    b(i,:) = mapnonzero(E(i),f,ktrial);
                end
                if isvector(b), b = b(:)'; end
                return
            end
            % singleton E
            b = false(1,length(ktrial));
            % checks
            [fname siz offset] = mapfile(E,f);
            if isempty(fname), return, end
            % read file
            fid = fopen(fname,'r+');
            for i=1:length(ktrial)
                fseek(fid,offset + (ktrial-1)*prod(siz(1:2))*4,'bof');
                b(i) = logical(fread(fid,1,'single'));
            end
            fclose(fid);
        end
        function x = get.data(E)
            x = mapread(E,'data');
        end
        function x = getdata(E,varargin)
            % function x = getdata(E[,whichtrials][,'avg|'sum'][,classname])
            x = mapread(E,'data',varargin{:});
        end
        function x = get.V(E)
            x = mapread(E,'V');
        end
        function x = getV(E,varargin)
            x = mapread(E,'V',varargin{:});
        end
        function x = get.J(E)
            x = mapread(E,'J');
        end
        function x = getJ(E,varargin)
            x = mapread(E,'J',varargin{:});
        end
        function x = get.G(E)
            x = mapread(E,'G');
        end
        function x = get.R(E)
            x = mapread(E,'R');
        end
        function x = getG(E,varargin)
            x = mapread(E,'G',varargin{:});
        end
        function x = getR(E,varargin)
            x = mapread(E,'R',varargin{:});
        end
        function setdata(E,ktrial,x)
            mapwrite(E,'data',ktrial,x)
        end
        function setflow(E,ktrial,varargin)
            % for example, setflow(E,ktrial,V,J)
            % checks
            F = fieldnames(E.savepar);
            if length(varargin)~=length(F)
                error('wrong number of estimation results')
            end
            for i=1:length(F)
                mapwrite(E,F{i},ktrial,varargin{i})
            end
        end
        function b = hasdata(E,ktrial)
            b = mapnonzero(E,'data',ktrial);
        end
        function b = hasresult(E,ktrial)
            % multiple E?
            if length(E)>1
                b = false(length(E),length(ktrial));
                for i=1:length(E)
                    b(i,:) = hasresult(E(i),ktrial);
                end
                if isvector(b), b = b(:)'; end
                return
            end
            % singleton E
            F = fieldnames(E.savepar);
            b = mapnonzero(E,F{1},ktrial);
        end
        function cleanfiles(E)
            % function cleanfiles(E)
            %---
            % remove files linked to edge E with obsolete keys
            nE = length(E);
            all = cell(1,nE);
            keep = cell(1,nE);
            for k=1:length(numel(E))
                e = E(k);
                d = dir([e.savedir 'edge' e.id_edge '*']);
                all{k} = {d.name};
                fdata = mapfile(e,'data','create'); [dum fdata] = fileparts(fdata);
                F = fieldnames(e.savepar);
                fflow = mapfile(e,F{1},'create'); [dum fflow] = fileparts(fflow);
                keep{k} = {fdata fflow};
            end
            all = [all{:}];
            keep = [keep{:}];
            rem = setdiff(all,keep);
            if isempty(rem), disp('no need to delete file'), return, end
            confirmation = questdlg([{'Are you sure to delete files?'},rem], ...
                'Confirm file deletion','Yes','No','No');
            if strcmp(confirmation,'Yes')
                swd = pwd; cd(e.savedir)
                delete(rem{:})
                cd(swd)
            end
        end
    end
    
    % other routines
    methods
        function[tubea tubeb] = gettube(E)
            dsnake = gradient(E.snake);
            dnorm  = sqrt(sum(dsnake.^2));
            snaken = [dsnake(2,:)./dnorm; -dsnake(1,:)./dnorm];
            hdiam = repmat(E.halfd,2,1);
            tubea = E.snake + snaken.*hdiam;
            tubeb = E.snake - snaken.*hdiam;
        end
        function y = getdataf(E,k)
            y = getdata(E,k,'double');
            for i=1:length(E.filterpar)
                eval(E.filterpar{i});
            end
        end
        function b = checksavedir(E)
            b = exist(E.savedir,'dir');
            if ~b, disp('Please set save directory'), end
        end
        function releasemap(E)
            E.datamap = [];
            E.flowmap = [];
        end
    end
    
    
end



