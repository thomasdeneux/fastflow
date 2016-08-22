classdef fast_data < hgsetget
    
    properties (SetAccess='protected')
        % fastflow object
        S
        % temporary data
        tmpdata = struct;
    end
    properties (Access='protected', Transient)
        tmpdatahidden
    end
    properties
        % user
        active = true;
        note = '';
        user
    end
    properties (SetAccess='protected')
        % file info
        filebase
        algo
        id
        % size
        np
        nt
        nexp
    end
    properties (Dependent, SetAccess='private')
        savedir
        files
    end
    
    % Constructor, copy
    methods
        function D = fast_data(S)
            % function D = fast_data(S)
            if nargin==0, return, end
            
            % set properties
            D.S = S;
            % note that np, nt, id, filename and algo remain to be set!
        end
    end
    
    % Data: file access
    methods (Access='protected')
        function [fname siz] = mapfile(D,docreate,linesflag)
            % function fname = mapfile(D,docreate[,'_suffix'])
            % function fname = mapfile(D,docreate[,{'_suffix',customsize}])
            siz = [D.np D.nt D.nexp];
            if nargin<2, docreate = false; end
            if nargin<3
                linesflag = ''; 
            elseif iscell(linesflag)
                specialsize = linesflag{2};
                siz(1:length(specialsize)) = specialsize;
                linesflag = linesflag{1};
            end
            fil = [D.filebase '_' D.algo];   
            if ~isempty(linesflag), fil = [fil '_' linesflag]; end
            fname = [D.savedir fil];
            %fprintf('access file %s\n',fname)
            if ~exist(fname,'file')
                if docreate
                    % create file if necessary
                    disp(['create file ' fil])
                    fid = fopen(fname,'w');
                    fwrite(fid,zeros(siz,'single'),'single');
                    fclose(fid);
                else
                    fname = '';
                end
            end
        end
        function x = mapread(D,ktrial,linesflag)
            % function x = mapread(D[,ktrial[,linesflag]])
            if nargin<3, linesflag = ''; end
            % file name and size
            [fname siz] = mapfile(D,false,linesflag);
            if any(ktrial<1 | ktrial>siz(3)), x=[]; return, end
            % Read file (in a block or trial by trial)
            if isempty(fname), x=[]; return, end
            fid = fopen(fname,'r');
            if isempty(ktrial) % not really faster to read a single block
                x = reshape(fread(fid,prod(siz),'single=>single'),siz);
            else
                x = zeros(siz(1),siz(2),length(ktrial),'single');
                for i=1:length(ktrial)
                    fseek(fid,(ktrial(i)-1)*prod(siz(1:2))*4,'bof');
                    x(:,:,i) = fread(fid,siz(1:2),'single=>single');
                end
            end
            fclose(fid);
        end
        function mapwrite(D,ktrial,x,linesflag)
            % function mapwrite(D,ktrial,x[,linesflag])
            if nargin<4, linesflag = ''; end    
            if ~isscalar(ktrial) || ktrial<1 || ktrial>D.nexp || mod(ktrial,1)
                error('wrong trial number')
            end
            % file name and size
            [fname siz] = mapfile(D,true,linesflag);
            if any(size(x)~=siz(1:2)), error('wrong size'), end
            % write in file
            fid = fopen(fname,'r+');
            fseek(fid,(ktrial-1)*prod(siz(1:2))*4,'bof');
            fwrite(fid,x,'single');
            fclose(fid);
            % mark that this trial has result now
            maphascontent(D,[],false,linesflag) % this is to create D.tmpdata.hasdata if necessary
            if iscell(linesflag)
                fieldname = ['hasdata_' linesflag{1}]; 
            elseif ~isempty(linesflag) 
                fieldname = ['hasdata_' linesflag];
            else
                fieldname = 'hasdata';
            end
            D.tmpdata.(fieldname)(ktrial) = true;
        end
        function settmpdata(D,ktrial,dataflag,x)
            % function settmpdata(D,ktrial,dataflag,x)
            savename = [dataflag '_current'];
            D.tmpdata.(savename).ktrial = ktrial;
            D.tmpdata.(savename).data   = x;
        end
        function b = maphascontent(D,ktrial,doforce,linesflag)
            % function b = maphascontent(D,ktrials[,doforce[,linesflag]])
            %---
            % if no output, only precomputes D.tmpdata.hasdata
            if nargin<3, doforce=false; end
            if nargin<4, linesflag=''; end
            if iscell(linesflag)
                fieldname = ['hasdata_' linesflag{1}]; 
            elseif ~isempty(linesflag) 
                fieldname = ['hasdata_' linesflag];
            else
                fieldname = 'hasdata';
            end
            
            % multiple D?
            if ~isscalar(D)
                b = false(length(D),length(ktrial));
                for i=1:length(D)
                    b(i,:) = maphascontent(D(i),ktrial,doforce,linesflag);
                end
                return
            end
            
            % compute and write in local storage tmpdata if not done yet
            if ~isfield(D.tmpdata,fieldname) || doforce
                % file name and size
                [fname siz] = mapfile(D,false,linesflag);
                if isempty(fname)
                    D.tmpdata.(fieldname) = false(1,D.nexp);
                else
                    % read file
                    fid = fopen(fname,'r+');
                    skip = (siz(1)*siz(2)-1) * 4; % read only one value from each trial
                    D.tmpdata.(fieldname) = logical(fread(fid,D.nexp,'single',skip));
                    fclose(fid);
                end
            end
            
            % output
            if nargout==1, b = D.tmpdata.(fieldname)(ktrial); end
        end
        function mapchgnexp(D,n)
            if n==D.nexp, return, end % nothing to do
            D.nexp = n;
            % increase size of associated files if required
            flist = getfiles(D,'full');
            if isempty(flist), return, end
            d = dir(flist{1});
            nold = d.bytes / (D.np*D.nt*4);
            if mod(nold,1), error programming, end
            if n>nold
                siz = [D.np D.nt n-nold];
                z = zeros(siz,'single');
                for i=1:length(flist)
                    fname = flist{i};
                    fid = fopen(fname,'a');
                    fwrite(fid,z,'single');
                end
            end
        end
    end
    methods
        function savedir = get.savedir(D)
            savedir = D.S.savedir;
        end
        function flist = getfiles(D,flag)
            % function flist = getfiles(D,'base|[full]')
            if nargin<2, flag='full'; end
            d = dir([D.savedir D.filebase '_' D.algo '*']);
            flist = {d.name};
            if strcmp(flag,'full') && ~isempty(flist)
                flist = strcat(D.savedir,flist);
            end
        end
        function flist = get.files(D)
            flist = getfiles(D,'base');
        end
    end
        
    % Results
    methods (Access='protected')
        function x = mapgetresult(D,resdisp)
            % function x = getresult(D,resdisp)
            % function x = getresult(D,ktrial)
            %---
            % Input:
            % - D           fast_data object
            % - resdisp     structure
            %
            % Output:
            % - x           array
            % 
            % The 'resdisp' structure can contain the following fields:
            %   * data      'volume', 'result', 'lines', 'section',
            %               'profile', 'sigma' [default is chosen according
            %               to the exact class of the fast_data object]
            %   * trial     trial number, ['average'], or 'timeavg' 
            %   * oddeven   ['all'], 'odd', or 'even'
            %   * res       ['all'], 'stim', 'rest', or 's/r'
            %   * points    scalar (smoothing), indices, ['image' or []],
            %               'all', 'middle', '2/5', '3/5', or '4/5' 
            %   * smooth    number
            %   * fzsub     integer number
            %
            % Output:
            % - x           3D array (space x time x trial) if several trials
            %               2D array (space x time-or-trial) if averaging
            %               column vector (time-or-trial) if also spatial averaging
            %               empty if no data available
            
            if nargin==1, help fastflow.getresult, x=[]; return, end
            
            % Input
            if ~isstruct(resdisp)
                ktrial = resdisp;
                resdisp = struct('trial',ktrial);
            end
            if ~isfield(resdisp,'data')
                switch class(D)
                    case 'fast_edge'
                        resdisp.data = 'volume';
                    case 'fast_flow'
                        resdisp.data = 'result';
                    case 'fast_section'
                        resdisp.data = 'width';
                end
            end
            if ~isfield(resdisp,'trial'), resdisp.trial = 'average'; end
            if ~isfield(resdisp,'res'), resdisp.res = 'stim/rest'; end
            if ~isfield(resdisp,'oddeven'), resdisp.oddeven = 'all'; end
                
            % Prepare
            % (volume scaling factor)
            if isfield(resdisp,'volumescale')
                volumescale = resdisp.volumescale; 
            else
                volumescale = 1;
            end
            % (frame zero subtraction)
            if isfield(resdisp,'fzsub')
                fzsub = resdisp.fzsub;
                if isfield(D.user,'tbin'), fzsub=ceil(fzsub/D.user.tbin); end
                fzsub = min(fzsub,D.nt);
            else
                fzsub = 0;
            end
            
            % Division between two conditions?
            if isequal(resdisp.trial,'average') && any(resdisp.res=='/')
                
                % Yes -> get the two conditions separately
                x = [];
                tokens = regexp(resdisp.res,'^(.*)/(.*)$','tokens'); tokens = tokens{1};
                resdisp.res = tokens{1}; x1 = getresultmain(D,resdisp); if isempty(x1), return, end
                resdisp.res = tokens{2}; x2 = getresultmain(D,resdisp); if isempty(x2), return, end
                
                % combination
                switch resdisp.data
                    case {'volume' 'section' 'flow'}
                        x = x1 ./ x2;
                        x = x * volumescale;
                    case 'sigma'
                        x = x1 ./ x2;
                    case 'result'
                        x = x1 - x2;
                    case {'volumef' 'lines'}
                        msgbox('what the hell are you trying to do!!??')
                    otherwise
                        error programming
                end
                if isempty(x), return, end
                
                % frame zero subtraction
                if fzsub
                    m = mean(x(:,1:fzsub),2); % x is 2D, m is 1D
                    switch resdisp.data
                        case {'volume' 'section' 'sigma' 'flow'}
                            x = fn_add(x,-m);
                        case 'result'
                            x = fn_add(x,-m);
                            
                            % normalize by baseline (not done yet because
                            % division subtraction was used instead of
                            % division)
                            m0 = abs(mean(mean((x1(:,1:fzsub))))+mean(mean((x2(:,1:fzsub)))))/2;
                            x = x / m0;
                    end
                end
                
                % normalize by peak (helps to plot together signals with
                % highly different dynamics)
                if strcmp(resdisp.data,'result')
                    peak = max(abs(mean(x,1)));
                    %x = x / peak;
                end
                
            else
                
                % No -> get the result straight
                x = getresultmain(D,resdisp);
                if isempty(x), return, end
                
                % frame zero subtraction
                if fzsub
                    m = mean(x(:,1:fzsub,:),2); % x is 2D or 3D
                    switch resdisp.data
                        case {'volume' 'section'}
                            x = fn_add(x,-m);
                            x = fn_mult(x,1./m);
                            x = x * volumescale;
                        case {'result' 'sigma' 'flow'}
                            x = fn_add(x,-m) / mean(m(:));
                        case {'volumef' 'lines'}
                            msgbox('what the hell are you trying to do!!??')
                            x = [];
                            return
                        otherwise
                            error programming
                    end
                end
                
            end
            
            % Column vector
            if isvector(x), x = x(:); end
            
        end
    end
    methods (Access='private')
        function x = getresultmain(D,resdisp) 
            % returns a 2D array (space x time-or-trial) 
            % or 3D array (space x time x trial)
            x = []; 
                      
            % Check on data type
            correctdata = fn_switch(class(D), ...
                'fast_edge',    {'volume' 'flow'}, ...
                'fast_section', {'section', 'profile', 'sigma'}, ...
                'fast_flow',    {'result', 'volumef', 'lines'});
            if ~isfield(resdisp,'data')
                resdisp.data = correctdata{1};
            else
                if ~ismember(resdisp.data,correctdata), error('data flag'), end
            end
            % (special: flow is computed from velocity and section width)
            if strcmp(resdisp.data,'flow')
                if isempty(D.flow) || isempty(D.section), return, end
                resdisp.data = 'result'; x1 = getresultmain(D.flow(1),resdisp);    if isempty(x1), return, end
                resdisp.data = 'sigma';  x2 = getresultmain(D.section(1),resdisp); if isempty(x2), return, end
                x = pi * x1.*((x2/2).^2);
                if size(x,1)>1, error('spatial averaging is mandatory for flow computation'), end
                return
            end
            % (special case for fast_section / sigma)
            if strcmp(resdisp.data,'sigma')
                % get the first line of the 'profile' result
                resdisp.data = 'profile';
                resdisp.points = 1; 
            end
                
            
            % Trials; use local storage or read directly the data
            if isnumeric(resdisp.trial)
                ktrial = resdisp.trial;
                if ~isfield(resdisp,'data')
                    resdisp.data = fn_switch(isa(D,'fast_edge'),'volume','result');
                end
                switch resdisp.data
                    case 'volumef'
                        x = getdataf(D,ktrial); % note that getdataf uses local store for data
                    otherwise
                        if isscalar(ktrial) 
                            % local storage
                            x = getlocalstore(D,resdisp.data,'trial',ktrial); 
                        else
                            % no local storage
                            linesflag = fn_switch(resdisp.data, ...
                                'lines','lines', ...
                                'profile',{'width',5}, ...
                                '');
                            x = mapread(D,ktrial,linesflag);
                        end
                end
            elseif strcmp(resdisp.trial,'timeavg')
                ktrial = []; % all trials
                switch resdisp.data
                    case 'volumef'
                        % who would ask for such thing?
                        x = getdataf(D,ktrial);
                        x = permute(mean(x,2),[1 3 2]);
                    otherwise
                        x = getlocalstore(D,resdisp.data,'timeavg');
                end
            else % 'average'
                % which conditions?
                if isempty(D.S.cond) || isscalar(D.S.cond)
                    icond = 0;
                elseif strcmp(resdisp.res,'all')
                    icond = unique(D.S.cond);
                elseif strcmp(resdisp.res,'rest')
                    icond = 0;
                elseif strcmp(resdisp.res,'stim')
                    icond = setdiff(unique(D.S.cond),0);
                else
                    token = regexp(resdisp.res,'^stim(\d+)$','tokens');
                    if isempty(token), error('wrong condition specification'), end
                    icond = str2double(token{1}{1});
                end
                ncond = length(icond);
                if ncond==0, return, end % can happen if no stimulated condition for example
                % which of odd and even?
                if strcmp(resdisp.oddeven,'all')
                    oddeven = {'odd' 'even'};
                else
                    oddeven = {resdisp.oddeven};
                end
                noddeven = length(oddeven);
                % use local storage?
                if strcmp(resdisp.data,'volumef')
                    % who would ask for such thing?
                    if isempty(D.S.cond) || isscalar(D.S.cond)
                        ktrial=1:D.nexp; 
                    else
                        ktrial = find(D.S.trials & ismember(D.S.cond,icond));
                    end
                    switch resdisp.oddeven
                        case 'odd'
                            ktrial = ktrial(1:2:end);
                        case 'even'
                            ktrial = ktrial(2:2:end);
                    end
                    x = getdataf(D,ktrial);
                    b = (x(1,1,:)~=0);
                    x = mean(x(:,:,b),3);
                else
                    % use local storage, need to gather several results!
                    ntrial = 0;
                    x = 0;
                    for i=1:ncond
                        for j=1:noddeven
                            [xij nij] = getlocalstore(D,resdisp.data,'average',icond(i),oddeven{j});
                            if nij
                                ntrial = ntrial+nij;
                                x = x + xij*nij; 
                            end
                        end
                    end
                    if ntrial
                        x = x/ntrial;
                    else
                        x = []; 
                        return
                    end
                end
            end
            if isempty(x), return, end
            
            % Finish
            x = processresult(D,x,resdisp);
        end
        function [x ntrial] = getlocalstore(D,dataflag,avgflag,varargin)
            % function [x ntrial] = getlocalstore(D,dataflag,'trial',ktrial)
            % function [x ntrial] = getlocalstore(D,dataflag,'timeavg')
            % function [x ntrial] = getlocalstore(D,dataflag,'average',icond,oddeven)
            %---
            % attempt to get result from the local storage; if not
            % successfull, reads from the main result files, computes, and
            % stores in the local storage
            % syntax 1 -> single trial, result or lines
            % syntax 2 -> one point per trial
            % syntax 3 -> average over trials according to condition and odd/even
            
            % Default result
            x = []; ntrial = 0;
            
            % Input and saving name
            switch avgflag
                case 'trial'
                    if nargin~=4, error('number of arguments'), end
                    ktrial = varargin{1};
                    if ~isscalar(ktrial), error('only one trial for local store'), end
                    savename = [dataflag '_current'];
                case 'timeavg'
                     if nargin~=3, error('number of arguments'), end                   
                    savename = [dataflag '_timeavg'];
                case 'average'
                    if nargin~=5, error('number of arguments'), end
                    [icond oddeven] = deal(varargin{:});
                    savename = [dataflag '_cond' num2str(icond) oddeven];
                otherwise
                    error argument
            end
            
            % Try to read result from local storage tmpdata 
            if isempty(D.tmpdata), D.tmpdata=struct; end
            if isfield(D.tmpdata,savename) 
                if strcmp(avgflag,'trial')
                    if D.tmpdata.(savename).ktrial==ktrial
                         x = D.tmpdata.(savename).data;
                         return
                    end
                else
                    x = D.tmpdata.(savename).data;
                    ntrial = D.tmpdata.(savename).ntrial;
                    return
                end                    
            end
            
            % Trials
            switch avgflag
                case 'trial'
                    % ktrial already set
                case 'timeavg'
                    ktrial = 1:D.nexp; % all trials
                case 'average'
                    [icond oddeven] = deal(varargin{:});
                    if isempty(D.S.cond) || isscalar(D.S.cond)
                        ktrial = 1:D.nexp;
                    else
                        ktrial = find(D.S.trials & ismember(D.S.cond,icond));
                    end
                    switch oddeven
                        case 'odd'
                            ktrial = ktrial(1:2:end);
                        case 'even'
                            ktrial = ktrial(2:2:end);
                    end
                    if isempty(ktrial), return, end
            end
            
            % Read file 
            linesflag=fn_switch(dataflag, ...
                'lines','lines', ...
                'profile',{'width',5}, ...
                '');
            x = mapread(D,ktrial,linesflag);
            if isempty(x), return, end
            
            % Update trials (some might be missing because estimation has
            % not been performed yet)
            b = (x(1,1,:)~=0);
            if ~any(b), x = []; return, end
            ktrial = ktrial(b);
            ntrial = length(ktrial);
            
            % Averaging over trials or time
            switch avgflag
                case ''
                    % nothing to do
                case 'timeavg'
                    x = permute(mean(x,2),[1 3 2]);
                case 'average'
                    x = x(:,:,b); 
                    x = mean(x,3);
            end
            
            % Store in local storage tmpdata
            if strcmp(avgflag,'trial')
                s = struct('ktrial',ktrial,'data',x);
            else
                s = struct('ntrial',ntrial,'data',x);
            end
            D.tmpdata.(savename) = s;           
        end
        function x = processresult(D,x,resdisp)
            % Spatial averaging
            if isfield(resdisp,'points')
                points = resdisp.points;
            else 
                points = [];
            end
            if ischar(points)
                switch points
                    case 'image'
                        points = [];
                    case 'all'
                        points = 1:D.np;
                    case 'middle'
                        points = floor(1+D.np/3):ceil(2*D.np/3);
                    case {'2/5' '3/5' '4/5'}
                        k5 = str2double(points(1));
                        points = floor(1+(k5-1)*D.np/5):ceil(k5*D.np/5);
                    otherwise
                        error('wrong point flag ''%s''',points)
                end
            end
            if isscalar(points) && points<0
                % spatial smoothing
                if isfield(D.user,'xbin'), points=-points/D.user.xbin; end
                x = filty(x,-points,'lm');
            elseif ~isempty(points)
                % change in syntax from previous version
                if isscalar(points) && (points~=1 || ~strcmp(class(D),'fast_section'))
                    error('change in syntax to get spatially smoothed result, please use a negative value for resdisp.points')
                end
                % spatial averaging
                x = mean(x(points,:,:),1);
            end
        
            % Smoothing
            if isfield(resdisp,'smooth') && resdisp.smooth
                sigma = resdisp.smooth;
                if isfield(D.user,'tbin'), sigma=sigma/D.user.tbin; end
                x = filtx(x,sigma,'lm');
            end
            
            %             % Scaling if required
            %             if strcmp(resdisp.data,'result') && isfield(D.user,'resultscale')
            %                 x = x*D.user.resultscale;
            %             end
        end
    end
    
    % Misc
    methods
        function chgfield(D,fieldname,value)
            if isempty(D), return, end
            [D.(fieldname)] = deal(value);
        end
        function hidetmpdata(D,flag)
            % function hidetmpdata(D,'remove|hide|show')
            if isempty(D), return, end
            switch flag
                case 'remove'
                    [D.tmpdata] = deal(struct);
                    [D.tmpdatahidden] = deal([]);
                case 'hide'
                    [D.tmpdatahidden] = deal(D.tmpdata);
                    [D.tmpdata] = deal(struct);
                case 'show'
                    [D.tmpdata] = deal(D.tmpdatahidden);
                    [D.tmpdatahidden] = deal([]);
            end
        end
        function copyfiles(D,D1,trialsoff,doflip)
            % function copyfiles(D,D1[,trialsoff[,doflip]])
            %---
            % copy in D the data/results from D1
            % set trialsoff to not copy from specific trials
            % set doflip to true if the vessel has been fliped upside-down
            if nargin<3, trialsoff=[]; end             
            if nargin<4, doflip=false; end
            if any([D.np D.nt D.nexp]~=[D1.np D1.nt D1.nexp]), error('size mismatch'), end
            if ~strcmp(D.method,D1.method), error('method mismatch'), end
            filesold = D1.files;
            for k=1:length(filesold)
                fold = filesold{k};
                fnew = strrep(fold,[D1.filebase '_' D1.algo],[D.filebase '_' D.algo]);
                disp([fold ' -> ' fnew])
                if doflip || any(trialsoff)
                    % more difficult
                    fid = fopen([D1.savedir fold],'r');
                    x = fread(fid,[D.np D.nt*D.nexp],'single=>single');
                    if doflip, 
                        x = flipud(x);
                        if ~any(findstr(fold,'volume')) && ~any(findstr(fold,'lines'))
                            % sign inversion for velocity result
                            x = -x;
                        end
                    end
                    if any(trialsoff)
                        x = reshape(x,[D.np D.nt D.nexp]);
                        x(:,:,trialsoff) = 0;
                    end                       
                    fclose(fid);
                    fid = fopen([D.savedir fnew],'w');
                    fwrite(fid,x,'single');
                    fclose(fid);
                    % delete(fold)
                else
                    copyfile([D1.savedir fold],[D.savedir fnew])
                    % movefile(fold,fnew)
                end
            end
        end
    end
    
    
end
