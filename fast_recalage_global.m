function X = fast_recalage_global(varargin)
% function X = fast_recalage_global(CS[,method][,x0][,'noresample'])
% function fast_recalage_global([method,]X)
%---
%
% Global input and output:
% - Y       3D array - movie to be coregistered and/or resampled
% 
% Input:
% - CS      reference image
% - X       parameters of the coregistration for the whole movie
% - x0      initialization parameters (for the whole movie or the first
%           frame, depending on the method used)
% - method  available methods are: 'isometry', 'nonlinear'
% - 'noresample'    use this flag for no resample
% 
% Output:
% - X       parameters of the coregistration for the whole movie

if nargin==0, help fast_recalage_global , return, end 

global Y 

if isempty(Y)
    error('global variable is not set')
end

[ni nj nt] = size(Y);


% Input
a = varargin{1};
doestimate = (isnumeric(a) && isequal(size(a),[ni nj]));
doresample = true;
method = 'isometry';
if doestimate
    CS = a;
    x = [];
    for i=2:length(varargin)
        a = varargin{i};
        if ischar(a)
            if strcmp(a,'noresample')
                doresample = false;
            else
                method = a;
            end
        else
            x = a;
        end
    end
else
    for i=1:length(varargin)
        a = varargin{i};
        if ischar(a)
            method = a;
        else
            X = a;
        end
    end
end

% Estimation
if doestimate
    % handle class
    switch class(Y)
        case 'single'
            CS = single(CS);
        case 'double'
            CS = double(CS);
    end
    
    switch method
        case {'isometry' 'isometry(frame)'}
            fn_progress('registering:',nt)
            X = zeros(nt,3);
            par = fn_register('par');
            par.ref = CS;
            for i=1:nt
                fn_progress(i)
                if ~isempty(x), par.shift0 = x(end-1:end); end
                x = fn_register(Y(:,:,i),par);                
                X(i,2:3) = x'; % no more rotation!!!
            end
        case  'isometry(trial)'
            disp('registering average frame')
            % high pass!
            CS2 = filt2(CS,15,'h');
            mY2 = filt2(mean(Y,3),15,'h');
            par = fn_register('par');
            par.ref = CS2;
            if ~isempty(x), par.shift0 = x(end-1:end); end
            X = zeros(nt,3);
            X(:,2:3) = fn_register(mY2,par)';
        case {'nonlinear' 'nonlinear(trial)'}
            X = fast_register_nonlinear(CS,mean(Y,3),x);
        otherwise
            error('unknown method ''%s''',method)
    end
end

% Resampling
if doresample
    switch method
        case {'isometry' 'isometry(frame)'}            
            fn_progress('resampling',nt)
            for i=1:nt
                fn_progress(i)
                Y(:,:,i) = fast_register(Y(:,:,i),X(i,:));
            end
        case  'isometry(trial)'
            fn_progress('resampling',nt)
            for i=1:nt
                fn_progress(i)
                Y(:,:,i) = fast_register(Y(:,:,i),X);
            end
        case {'nonlinear' 'nonlinear(trial)'}
            xplus = fast_register_nonlinear(Y(:,:,1),X,'parameters');
            fn_progress('resampling',nt)
            for i=1:nt
                if ~mod(i,10), fn_progress(i), end
                Y(:,:,i) = fast_register_nonlinear(Y(:,:,i),xplus);
            end            
        otherwise
            error('unknown method ''%s''',method)
    end
end

% Output
if nargout==0, clear X, end
