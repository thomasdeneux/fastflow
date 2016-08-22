classdef fast_flow < fast_data
    % function F = fast_flow(E,method,filterpar,estpar,id)
    
    properties (SetAccess='private')
        % parent fast_edge object
        edge
        method
        filterpar = {};
        binning = [1 1];
        estpar
    end
    properties (Dependent, SetAccess='private', Transient)
        result
        lines
    end
    
    % Constructor
    methods
        function F = fast_flow(E,method,filterpar,estpar,id)
            % function F = fast_flow(E,method,filterpar,estpar,id)
            if nargin==0, arg={}; else arg={E.S}; end
            F = F@fast_data(arg{:});
            if nargin==0, return, end
            if nargin<5, id = fast_flow.getflowid(method,filterpar,estpar); end
            F.edge = E;
            % spatial/temporal size can change! apply filter to check what
            % is the new size
            y = zeros(E.np,E.nt);
            for kfilt=1:length(filterpar)
                str = filterpar{kfilt};
                eval(str);
                if findstr(str,'fn_bin')
                    token = regexp(str,'fn_bin\(y,\[(\d+) +(\d+)\]\)','tokens');
                    if isempty(token) || isfield(F.user,'xbin')
                        error('binning detected, but could not interpret it')
                    end
                    F.user.xbin = str2double(token{1}{1});
                    F.user.tbin = str2double(token{1}{2});
                    F.user.resultscale = F.user.xbin/F.user.tbin;
                end
            end
            [F.np F.nt] = size(y);
            F.nexp = E.nexp;
            F.method = method;
            F.filterpar = filterpar;
            F.estpar = estpar;
            F.filebase = E.filebase;
            F.id = id;    
            F.algo = [F.method F.id];
        end
    end
    
    % Results: file access
    methods
        function x = get.result(D)
            x = mapread(D,[]);
        end
        function x = get.lines(D)
            x = mapread(D,[],'lines');
        end
        function x = getresult(F,a)
            % function x = getresult(F,ktrials)
            % function x = getresult(F,resdisp)
            if isstruct(a)
                resdisp = a;
            else
                resdisp = struct('data','result','trial',a);
            end
            if ~isscalar(F)
                x = cell(size(F));
                for k=1:numel(F)
                    x{k} = mapgetresult(F(k),a);
                end
            else
                x = mapgetresult(F,resdisp);
            end
        end
        function x = getlines(D,a)
            % function x = getlines(E,ktrials)
            % function x = getlines(E,resdisp)
            if isstruct(a)
                resdisp = a;
                resdisp.data = 'lines';
            else
                resdisp = struct('data','lines','trial',a);
            end
            x = mapgetresult(D,resdisp);
        end
        function b = hasresult(F,ktrial,doforce)
            % function b = hasresult(F,ktrials[,doforce])
            if nargin<2, ktrial=1:F.nexp; end
            if nargin<3, doforce=false; end
            b = maphascontent(F,ktrial,doforce);
        end
        function setresult(F,ktrial,x,y)
            mapwrite(F,ktrial,x)
            if nargin==4, mapwrite(F,ktrial,y,'lines'), end
        end        
        function setresulttmp(F,ktrial,x,y)
            settmpdata(F,ktrial,'result',x);
            if nargin==4, settmpdata(F,ktrial,'lines',y), end
        end        
    end
    
    % Estimation: data filtering
    methods
        function z = getdataf(F,ktrial)
            % function z = getdataf(F,ktrial)
            %---
            % note that this function does not admit a syntax with the
            % resdisp structure, as it does not apply any averaging or
            % other postprocessing to the data
            x = getdata(F.edge,ktrial);
            x = double(x);
            if isempty(x), z=[]; return, end
            n = size(x,3);
            z = zeros([F.np F.nt n],'single');
            for i=1:n
                y = x(:,:,i);
                for kfilt=1:length(F.filterpar)
                    eval(F.filterpar{kfilt});
                end
                z(:,:,i) = y;
            end
        end
    end
    
    % flow ID
    methods (Static)
        function [idflow s] = getflowid(method,filterpar,estpar)
            s = struct('METHOD',method, ...
                'FILTER',{filterpar(:)},'ESTPAR',estpar);
            idflow = fn_hash(s,4);
        end
    end
    
    % Misc
    methods
        function access(F) %#ok<MANU>
            keyboard
        end
        function eval(F,command) %#ok<MANU>
            eval(command)
        end
    end
    
end