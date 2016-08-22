classdef fast_section < fast_data
    % function U = fast_section(E,sectionpar,id)
    
    properties (SetAccess='private')
        % parent edge
        edge
        % section parameters
        sectionpar
        dx
    end
    properties (Dependent, SetAccess='private', Transient)
        % section data
        sectiondata
        % result
        gaussianprofile
        sigma
    end
    
    % Constructor
    methods
        function U = fast_section(E,sectionpar,id)
            % function F = fast_section(E,sectionpar,estpar,id)
            if nargin==0, arg={}; else arg={E.S}; end
            U = U@fast_data(arg{:});
            if nargin==0, return, end
            if nargin<3, id = fast_section.getsectionid(sectionpar); end
            U.edge = E;
            U.S = E.S;
            x = fast_vesselsection(zeros(E.S.nx,E.S.ny),E,sectionpar); % test interpolation to determine number of points
            U.np = size(x,1);
            U.nt = E.nt;
            U.nexp = E.nexp;
            U.sectionpar = sectionpar;
            U.filebase = E.filebase;
            U.id = id;    
            U.algo = ['section' id];
        end
    end
    
    % Section data: file access
    methods
        function x = getsectiondata(U,a)
            % function x = getsectiondata(U,ktrials)
            % function x = getsectiondata(U,resdisp)
            if isstruct(a)
                resdisp = a;
            else
                resdisp = struct('data','section','trial',a);
            end
            x = mapgetresult(U,resdisp);
        end
        function setsectiondata(U,ktrial,x)
            % function setsectiondata(E,ktrial,x)
            %---
            % x is an np x nt array
            mapwrite(U,ktrial,x)
            if isfield(U.tmpdata,'current_section') && U.tmpdata.current_volume.ktrial==ktrial
                U.tmpdata = rmfield(U.tmpdata,'current_section');
            end 
        end
        function setsectiondatatmp(U,ktrial,x)
            % function setsectiondatatmp(U,ktrial,x)
            settmpdata(U,ktrial,'section',x)
        end
        function b = hassectiondata(U,ktrial,doforce)
            % function b = hassectiondata(U,ktrial[,doforce])
            if nargin<3, doforce=false; end
            b = maphascontent(U,ktrial,doforce);
        end
        function x = get.sectiondata(U)
            x = mapread(U,[]);
        end
    end
    
    % Result data: file access
    methods
        function setresult(U,ktrial,x)
            mapwrite(U,ktrial,x,{'width',[5 U.nt]})
        end
        function setresulttmp(U,ktrial,x)
            % function setresulttmp(U,ktrial,x)
            settmpdata(U,ktrial,'profile',x)
        end
        function b = hasresult(U,ktrial,doforce)
            % function b = hasresult(U,ktrial,doforce)
            if nargin<3, doforce=false; end
            b = maphascontent(U,ktrial,doforce,{'width',[5 U.nt]});
        end
        function x = getresult(U,a)
            % function x = getresult(U,ktrials)
            % function x = getresult(U,resdisp)
            if isstruct(a)
                resdisp = a;
            else
                resdisp = struct('data','sigma','trial',a);
            end
            if ~isscalar(U)
                x = cell(size(U));
                for k=1:numel(U)
                    x{k} = mapgetresult(U(k),a);
                end
            else
                x = mapgetresult(U,resdisp);
            end
        end
        function x = get.gaussianprofile(U)
            x = mapread(U,[],{'width',[5 U.nt]});
        end
        function x = get.sigma(U)
            x = U.profile(1,:,:);
        end
        function chgnexp(E,n) %#ok<MANU,INUSD>
            error('not implemented yet')
        end
    end
    
    % ID
    methods (Static)
        function idsection = getsectionid(sectionpar)
            idsection = fn_hash(sectionpar,4);
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