function res = fast_estimateflow(I,estpar,ha)
% function res = fast_estimateflow(I,estpar[,ha])
%---
% returns a structure with the results

persistent gaborpar gabornp gabornt gaborbase 
if isempty(gabornp), gabornp=0; gabornt=0; end

% Input
if nargin<3
    ha = {};
else
    ha = {ha};
end

% Radon
if findstr(estpar.method,'radon')
    par = fast_radon('par');
    estpar = fn_structmerge(par,estpar);
    res.R = fast_radon(I,estpar);
end

% Gabor
if findstr(estpar.method,'gabor')
    par = fast_gabor('par');
    estpar = fn_structmerge(par,estpar);

    % compute filter base if necessary
    tmp = fn_structmerge(fast_gabor('par'),estpar,'skip','type');
    [np nt] = size(I);
    if ~isequal(tmp,gaborpar) || np>gabornp
        gaborpar = tmp;
        gabornp = max(gabornp,np);
        gabornt = max(gabornt,nt);
        gaborbase = fast_gabor(gaborpar,gabornp,gabornt);
    end
    
    % estimate
    gabor = fast_gabor(I,gaborpar,gaborbase); % gaborpar and gaborbase defined before the trials loop
    switch estpar.useresult
        case 'anglesc'
            res.G = cotd(gabor.anglesc);
        otherwise
            error('unknown useresult flag ''%s'' to choose Gabor estimation result',estpar.useresult)
    end
end

% fast track
if findstr(estpar.method,'track')
    par = fast_track('par');
    estpar = fn_structmerge(par,estpar);
    
    % initialization
    flag = regexp(estpar.method,'gabor(.)track','tokens');
    if ~isempty(flag) && any(flag=='-+')
        error('usage is deprecated; use v0=''gabor'' to initialize with the mean of Gabor result')
    end
    if ischar(par.v0)
        switch lower(par.v0)
            case 'gabor'
                par.v0 = mean(res.G(:));
            case 'radon'
                par.v0 = mean(res.R(:));
        end
    end
    
    % estimation
    [res.J res.V] = fast_track(I,estpar,ha{:});

end


