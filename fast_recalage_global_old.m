function varargout = fast_recalage_global(varargin)
% function ISO = fast_recalage_global(CS[,iso0][,indices][,'noresample'])
% function fast_recalage_global(ISO)
%---
% achieves vessel registration (estimation of movement based on CS)
% and normalization (divide by mean)

if nargin==0, help fast_recalage_global , return, end 

global Y 

if isempty(Y)
    error('global variable is not set')
end

[nj ni nt] = size(Y);

a = varargin{1};

% Resampling only
if any(size(a)==3);
    ISO = a;
    if size(ISO,2)~=3, ISO=ISO'; end
    fn_progress('resampling',nt)
    for i=1:nt
        fn_progress(i)
        Y(:,:,i) = fast_register(Y(:,:,i),ISO(i,:));
    end
    varargout = {Y};
    return
end

% Coregister: input
CS = a;
iso = [0 0 0];
indices = [];
doresample = true;
for i=2:length(varargin)
    a = varargin{i};
    if ischar(a)
        if ~strcmp(a,'noresample'), error argument, end
        doresample = false;
    else
        if ~isvector(a), error argument, end
        if length(a)==3
            iso = a;
        else
            indices = a;
        end
    end
end
if isempty(indices)
    skipe = 20;
    bin = 3;
    [I,J] = meshgrid(1+skipe:bin:ni-skipe,1+skipe:bin:nj-skipe);
    indices = J(:) + nj*(I(:)-1);
end

% slight smoothing parameter
sigma = 2;

fn_progress('registering:',nt)
% register first frame to reference image using full window
fn_progress(1)
iso = fast_register(CS,Y(:,:,1),iso);
[nj ni] = size(CS);
M1 = fast_register(iso,nj,ni);
% register all frames to the first frame using smaller window (defined in
% 'indices') ; this way, even if the registering to reference image is
% problematic, frames remain registered together 
ISO = zeros(nt,3); iso = [0 0 0];
Y1 = filt2(Y(:,:,1),sigma);
for i=1:nt
    fn_progress(i)
    
    % register every frame to first frame
    iso = fast_register(Y1,filt2(Y(:,:,i),sigma),iso,indices);
    % deduce the registration of the frame to reference image
    M2 = fast_register(iso,nj,ni);
    ISO(i,:) = fast_register(M2*[M1;0 0 1],nj,ni);
    
    if doresample, Y(:,:,i) = fast_register(Y(:,:,i),ISO(i,:)); end
end
varargout = {ISO};

